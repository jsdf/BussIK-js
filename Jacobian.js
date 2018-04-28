// @flow
/*
*
* Inverse Kinematics software, with several solvers including
* Selectively Damped Least Squares Method
* Damped Least Squares Method
* Pure Pseudoinverse Method
* Jacobian Transpose Method
*
*
* Author: Samuel R. Buss, sbuss@ucsd.edu.
* Web page: http://www.math.ucsd.edu/~sbuss/ResearchWeb/ikmethods/index.html
*
*
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*
*
*/

import {VectorR3, VectorR3_Zero} from './LinearR3';
import VectorRn from './VectorRn';
import Node from './Node';
import Tree from './Tree';

type UpdateMode =
  | 'JACOB_Undefined'
  | 'JACOB_JacobianTranspose'
  | 'JACOB_PseudoInverse'
  | 'JACOB_DLS'
  | 'JACOB_SDLS';

// declare class Tree {
//   Init(): void;
//   Compute(): void;
//   GetRoot(): Node;
//   InsertRoot(root: Node): void;
//   InsertLeftChild(parent: Node, child: Node): void;

//   Compute(): void;
//   GetNumEffector(): number;
//   GetNumJoint(): number;
//   GetParent(n: Node): Node;
//   GetRoot(): Node;
//   GetSuccessor(n: Node): Node;
// }

type Ptr = {
  postincrementDeref(): number,
};

declare class MatrixRmn {
  SetSize(row: number, col: number): void;
  SetZero(): void;
  SetTriple(i: number, j: number, fill: VectorR3): void;
  MultiplyTranspose(v: VectorRn, result: VectorRn): void;
  Multiply(v: VectorRn, result: VectorRn): void;
}

function Square(x) {
  return x * x;
}

function assert(pass: boolean) {
  if (!pass) throw new Error('assertion failed');
}

const DegreesToRadians = Math.PI / 180;

function VectorRnDot(u: VectorRn, v: VectorRn) {
  assert(u.GetLength() == v.GetLength());
  let res = 0.0;
  const p: Ptr = u.GetPtr();
  const q: Ptr = v.GetPtr();
  for (let i = u.GetLength(); i > 0; i--) {
    res += p.postincrementDeref() * q.postincrementDeref();
  }
  return res;
}

// #ifdef _DYNAMIC
// const  BASEMAXDIST = 0.02;
// #else
const MAXDIST = 0.08; // optimal value for double Y shape : 0.08
// #endif
const DELTA = 0.4;
const LAMBDA = 2.0; // only for DLS. optimal : 0.24
const NEARZERO = 0.0000000001;

// Parameters for pseudoinverses
const PseudoInverseThresholdFactor = 0.01; // Threshold for treating eigenvalue as zero (fraction of largest eigenvalue)

// Parameters for damped least squares
// Optimal damping values have to be determined in an ad hoc manner  (Yuck!)
const DefaultDampingLambda = 0.6; // Optimal for the "Y" shape (any lower gives jitter)
//const DefaultDampingLambda = 1.1;      // Optimal for the DLS "double Y" shape (any lower gives jitter)
// const DefaultDampingLambda = 0.7;     // Optimal for the DLS "double Y" shape with distance clamping (lower gives jitter)

// Cap on max. value of changes in angles in single update step
const MaxAngleJtranspose = 30.0 * DegreesToRadians;
const MaxAnglePseudoinverse = 5.0 * DegreesToRadians;
const MaxAngleDLS = 45.0 * DegreesToRadians;
const MaxAngleSDLS = 45.0 * DegreesToRadians;
const BaseMaxTargetDist = 0.4;

export default class Jacobian {
  // public:
  // void ComputeJacobian(VectorR3 targets);
  //   void SetJendTrans(MatrixRmn J);
  //   void SetDeltaS(VectorRn S);

  // void CalcDeltaThetas();     // Use this only if the Current Mode has been set.
  // void ZeroDeltaThetas();
  // void CalcDeltaThetasTranspose();
  // void CalcDeltaThetasPseudoinverse();
  // void CalcDeltaThetasDLS();
  //   void CalcDeltaThetasDLS2(const VectorRn dVec);
  // void CalcDeltaThetasDLSwithSVD();
  // void CalcDeltaThetasSDLS();
  //   void CalcDeltaThetasDLSwithNullspace( const VectorRn desiredV);

  // void UpdateThetas();
  //   void UpdateThetaDot();
  // double UpdateErrorArray(VectorR3 targets);   // Returns sum of errors
  // void UpdatedSClampValue(VectorR3 targets);
  // void Reset();
  //
  // static void CompareErrors( const Jacobian& j1, const Jacobian& j2, double* weightedDist1, double* weightedDist2 );
  // static void CountErrors( const Jacobian& j1, const Jacobian& j2, int* numBetter1, int* numBetter2, int* numTies );

  ActiveJacobian(): MatrixRmn {
    return this.Jactive;
  }
  SetJendActive(): void {
    this.Jactive = this.Jend;
  } // The default setting is this.Jend.
  SetJtargetActive(): void {
    this.Jactive = this.Jtarget;
  }
  GetErrorArray() {
    return this.errorArray;
  }
  SetCurrentMode(mode: UpdateMode) {
    this.CurrentUpdateMode = mode;
  }
  GetCurrentMode(): UpdateMode {
    return this.CurrentUpdateMode;
  }
  SetDampingDLS(lambda: number) {
    this.DampingLambda = lambda;
    this.DampingLambdaSq = lambda * lambda;
  }

  GetNumRows() {
    return this.nRow;
  }
  GetNumCols() {
    return this.nCol;
  }

  // public:
  m_tree: Tree; // tree associated with this Jacobian matrix
  m_nEffector: number; // Number of end effectors
  nJoint: number; // Number of joints
  nRow: number; // Total number of rows the real J (= 3*number of end effectors for now)
  nCol: number; // Total number of columns in the real J (= number of joints for now)

  Jend: MatrixRmn; // Jacobian matrix based on end effector positions
  Jtarget: MatrixRmn; // Jacobian matrix based on target positions
  Jnorms: MatrixRmn; // Norms of 3-vectors in active Jacobian (SDLS only)

  U: MatrixRmn; // J = this.U * Diag(this.w) * this.V^T  (Singular Value Decomposition)
  w: VectorRn;
  V: MatrixRmn;

  CurrentUpdateMode: UpdateMode;

  dS: VectorRn; // delta s
  dT1: VectorRn; // delta t    --  these are delta S values clamped to smaller magnitude
  dSclamp: VectorRn; // Value to clamp magnitude of dT at.
  dTheta: VectorRn; // delta theta
  dPreTheta: VectorRn; // delta theta for single eigenvalue  (SDLS only)

  errorArray: VectorRn; // Distance of end effectors from target after updating

  // Parameters for damped least squares
  DampingLambda: number;
  DampingLambdaSq: number;
  //this.DampingLambdaSDLS:number;

  Jactive: MatrixRmn;

  constructor(tree: Tree) {
    this.m_tree = tree;
    this.m_nEffector = tree.GetNumEffector();
    this.nJoint = tree.GetNumJoint();
    this.nRow = 3 * this.m_nEffector; // Include only the linear part

    this.nCol = this.nJoint;

    this.Jend.SetSize(this.nRow, this.nCol); // The Jocobian matrix
    this.Jend.SetZero();
    this.Jtarget.SetSize(this.nRow, this.nCol); // The Jacobian matrix based on target positions
    this.Jtarget.SetZero();
    this.SetJendActive();

    this.U.SetSize(this.nRow, this.nRow); // The this.U matrix for SVD calculations
    this.w.SetLength(Math.min(this.nRow, this.nCol));
    this.V.SetSize(this.nCol, this.nCol); // The this.V matrix for SVD calculations

    this.dS.SetLength(this.nRow); // (Target positions) - (End effector positions)
    this.dTheta.SetLength(this.nCol); // Changes in joint angles
    this.dPreTheta.SetLength(this.nCol);

    // Used by Jacobian transpose method & DLS & SDLS
    this.dT1.SetLength(this.nRow); // Linearized change in end effector positions based on this.dTheta

    // Used by the Selectively Damped Least Squares Method
    //dT.SetLength(this.nRow);
    this.dSclamp.SetLength(this.m_nEffector);
    this.errorArray.SetLength(this.m_nEffector);
    this.Jnorms.SetSize(this.m_nEffector, this.nCol); // Holds the norms of the active J matrix

    this.Reset();
  }

  // constructor(bool useAngularJacobian,int nDof)
  // {

  //   this.m_tree = 0;
  //   this.m_nEffector = 1;

  //   if (useAngularJacobian)
  //   {
  //     this.nRow = 2 * 3 * this.m_nEffector; // Include both linear and angular part
  //   } else
  //   {
  //     this.nRow = 3 * this.m_nEffector; // Include only the linear part
  //   }

  //   this.nCol = nDof;

  //   this.Jend.SetSize(this.nRow, this.nCol);       // The Jocobian matrix
  //   this.Jend.SetZero();
  //   this.Jtarget.SetSize(this.nRow, this.nCol);      // The Jacobian matrix based on target positions
  //   this.Jtarget.SetZero();
  //   SetJendActive();

  //   this.U.SetSize(this.nRow, this.nRow);        // The this.U matrix for SVD calculations
  //   this.w .SetLength(Math.min(this.nRow, this.nCol));
  //   this.V.SetSize(this.nCol, this.nCol);        // The this.V matrix for SVD calculations

  //   this.dS.SetLength(this.nRow);     // (Target positions) - (End effector positions)
  //   this.dTheta.SetLength(this.nCol);   // Changes in joint angles
  //   this.dPreTheta.SetLength(this.nCol);

  //   // Used by Jacobian transpose method & DLS & SDLS
  //   this.dT1.SetLength(this.nRow);      // Linearized change in end effector positions based on this.dTheta

  //   // Used by the Selectively Damped Least Squares Method
  //   //dT.SetLength(this.nRow);
  //   this.dSclamp.SetLength(this.m_nEffector);
  //   this.errorArray.SetLength(this.m_nEffector);
  //   this.Jnorms.SetSize(this.m_nEffector, this.nCol);    // Holds the norms of the active J matrix

  //   Reset();
  // }

  Reset(): void {
    // Used by Damped Least Squares Method
    this.DampingLambda = DefaultDampingLambda;
    this.DampingLambdaSq = this.DampingLambda * this.DampingLambda;
    // this.DampingLambdaSDLS = 1.5*DefaultDampingLambda;

    this.dSclamp.Fill(Number.MAX_VALUE);
  }

  // Compute the deltaS vector, this.dS, (the error in end effector positions
  // Compute the J and K matrices (the Jacobians)
  ComputeJacobian(targets: Array<VectorR3>): void {
    if (this.m_tree == null) return;

    // Traverse tree to find all end effectors
    let temp: VectorR3;
    let n = this.m_tree.GetRoot();
    while (n) {
      if (n.IsEffector()) {
        let i: number = n.GetEffectorNum();
        const targetPos: VectorR3 = targets[i];

        // Compute the delta S value (differences from end effectors to target positions.
        temp = targetPos;
        temp.subtract(n.GetS());
        this.dS.SetTriple(i, temp);

        // Find all ancestors (they will usually all be joints)
        // Set the corresponding entries in the Jacobians J, K.
        let m = this.m_tree.GetParent(n);
        while (m) {
          let j: number = m.GetJointNum();
          assert(0 <= i && i < this.m_nEffector && 0 <= j && j < this.nJoint);
          if (m.IsFrozen()) {
            this.Jend.SetTriple(i, j, VectorR3_Zero);
            this.Jtarget.SetTriple(i, j, VectorR3_Zero);
          } else {
            temp = m.GetS(); // joint pos.
            temp.add(n.GetS()); // -(end effector pos. - joint pos.)
            temp.cross(m.GetW()); // cross product with joint rotation axis
            this.Jend.SetTriple(i, j, temp);
            temp = m.GetS(); // joint pos.
            temp.subtract(targetPos); // -(target pos. - joint pos.)
            temp.cross(m.GetW()); // cross product with joint rotation axis
            this.Jtarget.SetTriple(i, j, temp);
          }
          m = this.m_tree.GetParent(m);
        }
      }
      n = this.m_tree.GetSuccessor(n);
    }
  }

  // SetJendTrans(J: MatrixRmn): void {
  //   this.Jend.SetSize(J.GetNumRows(), J.GetNumColumns());
  //   this.Jend.LoadAsSubmatrix(J);
  // }

  // SetDeltaS(S: VectorRn): void {
  //   this.dS.Set(S);
  // }

  // The delta theta values have been computed in this.dTheta array
  // Apply the delta theta values to the joints
  // Nothing is done about joint limits for now.
  UpdateThetas(): void {
    // Traverse the tree to find all joints
    // Update the joint angles
    let n = this.m_tree.GetRoot();
    while (n) {
      if (n.IsJoint()) {
        let i: number = n.GetJointNum();
        n.AddToTheta(this.dTheta[i]);
      }
      n = this.m_tree.GetSuccessor(n);
    }

    // Update the positions and rotation axes of all joints/effectors
    this.m_tree.Compute();
  }

  UpdateThetaDot(): void {
    if (this.m_tree == null) return;

    // Traverse the tree to find all joints
    // Update the joint angles
    let n = this.m_tree.GetRoot();
    while (n) {
      if (n.IsJoint()) {
        let i: number = n.GetJointNum();
        n.UpdateTheta(this.dTheta[i]);
      }
      n = this.m_tree.GetSuccessor(n);
    }

    // Update the positions and rotation axes of all joints/effectors
    this.m_tree.Compute();
  }

  CalcDeltaThetas(): void {
    switch (this.CurrentUpdateMode) {
      case 'JACOB_Undefined':
        this.ZeroDeltaThetas();
        break;
      case 'JACOB_JacobianTranspose':
        this.CalcDeltaThetasTranspose();
        break;
      // case 'JACOB_PseudoInverse':
      //   this.CalcDeltaThetasPseudoinverse();
      //   break;
      // case 'JACOB_DLS':
      //   this.CalcDeltaThetasDLS();
      //   break;
      // case 'JACOB_SDLS':
      //   this.CalcDeltaThetasSDLS();
      //   break;
    }
  }

  ZeroDeltaThetas(): void {
    this.dTheta.SetZero();
  }

  // Find the delta theta values using inverse Jacobian method
  // Uses a greedy method to decide scaling factor
  CalcDeltaThetasTranspose(): void {
    const J: MatrixRmn = this.ActiveJacobian();

    J.MultiplyTranspose(this.dS, this.dTheta);

    // Scale back the this.dTheta values greedily
    J.Multiply(this.dTheta, this.dT1); // dT = J * this.dTheta
    const alpha: number = VectorRnDot(this.dS, this.dT1) / this.dT1.NormSq();
    assert(alpha > 0.0);
    // Also scale back to be have max angle change less than MaxAngleJtranspose
    const maxChange: number = this.dTheta.MaxAbs();
    const beta: number = MaxAngleJtranspose / maxChange;
    this.dTheta.multiply(Math.min(alpha, beta));
  }

  // CalcDeltaThetasPseudoinverse() :void
  // {
  //   const J:MatrixRmn  = this.ActiveJacobian();

  //   // Compute Singular Value Decomposition
  //   //  This an inefficient way to do Pseudoinverse, but it is convenient since we need SVD anyway

  //   J.ComputeSVD( this.U, this.w, this.V );

  //   // Next line for debugging only
  //     assert(J.DebugCheckSVD(this.U, this.w , this.V));

  //   // Calculate response vector this.dTheta that is the DLS solution.
  //   //  Delta target values are the this.dS values
  //   //  We multiply by Moore-Penrose pseudo-inverse of the J matrix
  //   const  pseudoInverseThreshold = PseudoInverseThresholdFactor*this.w.MaxAbs();

  //   const diagLength = this.w.GetLength();
  //   const wPtr: Ptr = this.w.GetPtr();
  //   this.dTheta.SetZero();
  //   for ( let i=0; i<diagLength; i++ ) {
  //     const dotProdCol = this.U.DotProductColumn( this.dS, i );    // Dot product with i-th column of this.U
  //     let alpha = wPtr.postincrementDeref();
  //     if ( fabs(alpha)>pseudoInverseThreshold ) {
  //       alpha = 1.0/alpha;
  //       MatrixRmn.AddArrayScale(this.V.GetNumRows(), this.V.GetColumnPtr(i), 1, this.dTheta.GetPtr(), 1, dotProdCol*alpha );
  //     }
  //   }

  //   // Scale back to not exceed maximum angle changes
  //   const maxChange = this.dTheta.MaxAbs();
  //   if ( maxChange>MaxAnglePseudoinverse ) {
  //     this.dTheta *= MaxAnglePseudoinverse/maxChange;
  //   }

  // }

  // CalcDeltaThetasDLSwithNullspace( desiredV:VectorRn):void
  // {
  //   const J:MatrixRmn = ActiveJacobian();

  //   MatrixRmn.MultiplyTranspose(J, J, this.U);    // this.U = J * (J^T)
  //   this.U.AddToDiagonal( this.DampingLambdaSq );

  //   // Use the next four lines instead of the succeeding two lines for the DLS method with clamped error vector e.
  //   // CalcdTClampedFromdS();
  //   // VectorRn dTextra(3*this.m_nEffector);
  //   // this.U.Solve( dT, &dTextra );
  //   // J.MultiplyTranspose( dTextra, this.dTheta );

  //   // Use these two lines for the traditional DLS method
  //   this.U.Solve( this.dS, this.dT1 );
  //   J.MultiplyTranspose( this.dT1, this.dTheta );

  //     // Compute JInv in damped least square form
  //     MatrixRmn UInv(this.U.GetNumRows(),this.U.GetNumColumns());
  //     this.U.ComputeInverse(UInv);
  //     assert(this.U.DebugCheckInverse(UInv));
  //     MatrixRmn JInv(J.GetNumColumns(), J.GetNumRows());
  //     MatrixRmn.TransposeMultiply(J, UInv, JInv);

  //     // Compute null space projection
  //     MatrixRmn JInvJ(J.GetNumColumns(), J.GetNumColumns());
  //     MatrixRmn.Multiply(JInv, J, JInvJ);
  //     MatrixRmn P(J.GetNumColumns(), J.GetNumColumns());
  //     P.SetIdentity();
  //     P -= JInvJ;

  //     // Compute null space velocity
  //     VectorRn nullV(J.GetNumColumns());
  //     P.Multiply(desiredV, nullV);

  //   // Compute residual
  //   VectorRn residual(J.GetNumRows());
  //   J.Multiply(nullV, residual);
  //   // TODO: Use residual to set the null space term coefficient adaptively.
  //   //printf("residual: %f\n", residual.Norm());

  //     // Add null space velocity
  //     this.dTheta += nullV;

  //   // Scale back to not exceed maximum angle changes
  //   double maxChange = this.dTheta.MaxAbs();
  //   if ( maxChange>MaxAngleDLS ) {
  //     this.dTheta *= MaxAngleDLS/maxChange;
  //   }
  // }

  // CalcDeltaThetasDLS():void
  // {
  //     const J:MatrixRmn = ActiveJacobian();

  //     MatrixRmn.MultiplyTranspose(J, J, this.U);    // this.U = J * (J^T)
  //     this.U.AddToDiagonal( this.DampingLambdaSq );

  //     // Use the next four lines instead of the succeeding two lines for the DLS method with clamped error vector e.
  //     // CalcdTClampedFromdS();
  //     // VectorRn dTextra(3*this.m_nEffector);
  //     // this.U.Solve( dT, &dTextra );
  //     // J.MultiplyTranspose( dTextra, this.dTheta );

  //     // Use these two lines for the traditional DLS method
  //     this.U.Solve( this.dS, &this.dT1 );
  //     J.MultiplyTranspose( this.dT1, this.dTheta );

  //     // Scale back to not exceed maximum angle changes
  //     double maxChange = this.dTheta.MaxAbs();
  //     if ( maxChange>MaxAngleDLS ) {
  //         this.dTheta *= MaxAngleDLS/maxChange;
  //     }
  // }

  // CalcDeltaThetasDLS2(const VectorRn dVec):void
  // {
  //     const J:MatrixRmn = ActiveJacobian();

  //     this.U.SetSize(J.GetNumColumns(), J.GetNumColumns());
  //     MatrixRmn.TransposeMultiply(J, J, this.U);
  //     this.U.AddToDiagonal( dVec );

  //     this.dT1.SetLength(J.GetNumColumns());
  //     J.MultiplyTranspose(this.dS, this.dT1);
  //     this.U.Solve(this.dT1, &this.dTheta);

  //     // Scale back to not exceed maximum angle changes
  //     double maxChange = this.dTheta.MaxAbs();
  //     if ( maxChange>MaxAngleDLS ) {
  //         this.dTheta *= MaxAngleDLS/maxChange;
  //     }
  // }

  // CalcDeltaThetasDLSwithSVD() :void
  // {
  //   const J:MatrixRmn = ActiveJacobian();

  //   // Compute Singular Value Decomposition
  //   //  This an inefficient way to do DLS, but it is convenient since we need SVD anyway

  //   J.ComputeSVD( this.U, this.w, this.V );

  //   // Next line for debugging only
  //     assert(J.DebugCheckSVD(this.U, this.w , this.V));

  //   // Calculate response vector this.dTheta that is the DLS solution.
  //   //  Delta target values are the this.dS values
  //   //  We multiply by DLS inverse of the J matrix
  //   long diagLength = this.w.GetLength();
  //   double* wPtr = this.w.GetPtr();
  //   this.dTheta.SetZero();
  //   for ( long i=0; i<diagLength; i++ ) {
  //     double dotProdCol = this.U.DotProductColumn( this.dS, i );    // Dot product with i-th column of this.U
  //     double alpha = *(wPtr++);
  //     alpha = alpha/(Square(alpha)+this.DampingLambdaSq);
  //     MatrixRmn.AddArrayScale(this.V.GetNumRows(), this.V.GetColumnPtr(i), 1, this.dTheta.GetPtr(), 1, dotProdCol*alpha );
  //   }

  //   // Scale back to not exceed maximum angle changes
  //   double maxChange = this.dTheta.MaxAbs();
  //   if ( maxChange>MaxAngleDLS ) {
  //     this.dTheta *= MaxAngleDLS/maxChange;
  //   }
  // }

  // CalcDeltaThetasSDLS() :void
  // {
  //   const J:MatrixRmn = ActiveJacobian();

  //   // Compute Singular Value Decomposition

  //   J.ComputeSVD( this.U, this.w, this.V );

  //   // Next line for debugging only
  //     assert(J.DebugCheckSVD(this.U, this.w , this.V));

  //   // Calculate response vector this.dTheta that is the SDLS solution.
  //   //  Delta target values are the this.dS values
  //   int this.nRows = J.GetNumRows();
  //   // TODO: Modify it to work with multiple end effectors.
  //   int numEndEffectors = 1;
  //   int this.nCols = J.GetNumColumns();
  //   this.dTheta.SetZero();

  //   // Calculate the norms of the 3-vectors in the Jacobian
  //   long i;
  //   const  *jx = J.GetPtr();
  //   double *jnx = this.Jnorms.GetPtr();
  //   for ( i=this.nCols*numEndEffectors; i>0; i-- ) {
  //     double accumSq = Square(*(jx++));
  //     accumSq += Square(*(jx++));
  //     accumSq += Square(*(jx++));
  //     *(jnx++) = sqrt(accumSq);
  //   }

  //   // Clamp the this.dS values
  //   CalcdTClampedFromdS();

  //   // Loop over each singular vector
  //   for ( i=0; i<this.nRows; i++ ) {

  //     double wiInv = this.w[i];
  //     if ( NearZero(wiInv,1.0e-10) ) {
  //       continue;
  //     }
  //     wiInv = 1.0/wiInv;

  //     double N = 0.0;           // N is the quasi-1-norm of the i-th column of this.U
  //     double alpha = 0.0;         // alpha is the dot product of dT and the i-th column of this.U

  //     const  *dTx = this.dT1.GetPtr();
  //     const  *ux = this.U.GetColumnPtr(i);
  //     long j;
  //     for ( j=numEndEffectors; j>0; j-- ) {
  //       double tmp;
  //       alpha += (*ux)*(*(dTx++));
  //       tmp = Square( *(ux++) );
  //       alpha += (*ux)*(*(dTx++));
  //       tmp += Square(*(ux++));
  //       alpha += (*ux)*(*(dTx++));
  //       tmp += Square(*(ux++));
  //       N += sqrt(tmp);
  //     }

  //     // M is the quasi-1-norm of the response to angles changing according to the i-th column of this.V
  //     //    Then is multiplied by the wiInv value.
  //     double M = 0.0;
  //     double *vx = this.V.GetColumnPtr(i);
  //     jnx = this.Jnorms.GetPtr();
  //     for ( j=this.nCols; j>0; j-- ) {
  //       double accum=0.0;
  //       for ( long k=numEndEffectors; k>0; k-- ) {
  //         accum += *(jnx++);
  //       }
  //       M += fabs((*(vx++)))*accum;
  //     }
  //     M *= fabs(wiInv);

  //     double gamma = MaxAngleSDLS;
  //     if ( N<M ) {
  //       gamma *= N/M;       // Scale back maximum permissable joint angle
  //     }

  //     // Calculate the this.dTheta from pure pseudoinverse considerations
  //     double scale = alpha*wiInv;     // This times i-th column of this.V is the psuedoinverse response
  //     this.dPreTheta.LoadScaled( this.V.GetColumnPtr(i), scale );
  //     // Now rescale the this.dTheta values.
  //     double max = this.dPreTheta.MaxAbs();
  //     double rescale = (gamma)/(gamma+max);
  //     this.dTheta.AddScaled(this.dPreTheta,rescale);
  //     /*if ( gamma<max) {
  //       this.dTheta.AddScaled( this.dPreTheta, gamma/max );
  //     }
  //     else {
  //       this.dTheta += this.dPreTheta;
  //     }*/
  //   }

  //   // Scale back to not exceed maximum angle changes
  //   double maxChange = this.dTheta.MaxAbs();
  //   if ( maxChange>MaxAngleSDLS ) {
  //     this.dTheta *= MaxAngleSDLS/(MaxAngleSDLS+maxChange);
  //     //this.dTheta *= MaxAngleSDLS/maxChange;
  //   }
  // }

  // CalcdTClampedFromdS() :void
  // {
  //   long len = this.dS.GetLength();
  //   long j = 0;
  //   for ( long i=0; i<len; i+=3, j++ ) {
  //     double normSq = Square(this.dS[i])+Square(this.dS[i+1])+Square(this.dS[i+2]);
  //     if ( normSq>Square(this.dSclamp[j]) ) {
  //       double factor = this.dSclamp[j]/sqrt(normSq);
  //       this.dT1[i] = this.dS[i]*factor;
  //       this.dT1[i+1] = this.dS[i+1]*factor;
  //       this.dT1[i+2] = this.dS[i+2]*factor;
  //     }
  //     else {
  //       this.dT1[i] = this.dS[i];
  //       this.dT1[i+1] = this.dS[i+1];
  //       this.dT1[i+2] = this.dS[i+2];
  //     }
  //   }
  // }

  // UpdateErrorArray(VectorR3 targets):double
  // {
  //   double totalError = 0.0;

  //   // Traverse tree to find all end effectors
  //   let temp:VectorR3;
  //   let n:Node = this.m_tree.GetRoot();
  //   while ( n ) {
  //     if ( n.IsEffector() ) {
  //       let i:number = n.GetEffectorNum();
  //       const targetPos:VectorR3 = targets[i];
  //       temp = targetPos;
  //       temp -= n.GetS();
  //       double err = temp.Norm();
  //       this.errorArray[i] = err;
  //       totalError += err;
  //     }
  //     n = this.m_tree.GetSuccessor( n );
  //   }
  //   return totalError;
  // }

  UpdatedSClampValue(targets: Array<VectorR3>): void {
    // Traverse tree to find all end effectors
    let temp: VectorR3;
    let n = this.m_tree.GetRoot();
    while (n) {
      if (n.IsEffector()) {
        let i: number = n.GetEffectorNum();
        const targetPos: VectorR3 = targets[i];

        // Compute the delta S value (differences from end effectors to target positions.
        // While we are at it, also update the clamping values in this.dSclamp;
        temp = targetPos;
        temp.subtract(n.GetS());
        let normSi = Math.sqrt(
          Square(this.dS[i]) + Square(this.dS[i + 1]) + Square(this.dS[i + 2])
        );
        let changedDist = temp.Norm() - normSi;
        if (changedDist > 0.0) {
          this.dSclamp[i] = BaseMaxTargetDist + changedDist;
        } else {
          this.dSclamp[i] = BaseMaxTargetDist;
        }
      }
      n = this.m_tree.GetSuccessor(n);
    }
  }

  // CompareErrors( const Jacobian& j1, const Jacobian& j2, double* weightedDist1, double* weightedDist2 ):void
  // {
  //   const e1:VectorRn = j1.this.errorArray;
  //   const e2:VectorRn = j2.this.errorArray;
  //   double ret1 = 0.0;
  //   double ret2 = 0.0;
  //   int len = e1.GetLength();
  //   for ( long i=0; i<len; i++ ) {
  //     double v1 = e1[i];
  //     double v2 = e2[i];
  //     if ( v1<v2 ) {
  //       ret1 += v1/v2;
  //       ret2 += 1.0;
  //     }
  //     else if ( v1 != 0.0 ) {
  //       ret1 += 1.0;
  //       ret2 += v2/v1;
  //     }
  //     else {
  //       ret1 += 0.0;
  //       ret2 += 0.0;
  //     }
  //   }
  //   *weightedDist1 = ret1;
  //   *weightedDist2 = ret2;
  // }

  // CountErrors( const Jacobian& j1, const Jacobian& j2, int* numBetter1, int* numBetter2, int* numTies ):void
  // {
  //   const e1:VectorRn = j1.this.errorArray;
  //   const e2:VectorRn = j2.this.errorArray;
  //   int b1=0, b2=0, tie=0;
  //   int len = e1.GetLength();
  //   for ( long i=0; i<len; i++ ) {
  //     double v1 = e1[i];
  //     double v2 = e2[i];
  //     if ( v1<v2 ) {
  //       b1++;
  //     }
  //     else if ( v2<v1 ) {
  //       b2++;
  //     }
  //     else {
  //       tie++;
  //     }
  //   }
  //   *numBetter1 = b1;
  //   *numBetter2 = b2;
  //   *numTies = tie;
  // }

  /* THIS VERSION IS NOT AS GOOD.  DO NOT USE!
CalcDeltaThetasSDLSrev2() :void 
{ 
  const J:MatrixRmn = ActiveJacobian();

  // Compute Singular Value Decomposition 

  J.ComputeSVD( this.U, this.w, this.V );
  
  // Next line for debugging only
    assert(J.DebugCheckSVD(this.U, this.w , this.V));

  // Calculate response vector this.dTheta that is the SDLS solution.
  //  Delta target values are the this.dS values
  int this.nRows = J.GetNumRows();
  int numEndEffectors = tree.GetNumEffector();   // Equals the number of rows of J divided by three
  int this.nCols = J.GetNumColumns();
  this.dTheta.SetZero();

  // Calculate the norms of the 3-vectors in the Jacobian
  long i;
  const  *jx = J.GetPtr();
  double *jnx = this.Jnorms.GetPtr();
  for ( i=this.nCols*numEndEffectors; i>0; i-- ) {
    double accumSq = Square(*(jx++));
    accumSq += Square(*(jx++));
    accumSq += Square(*(jx++));
    *(jnx++) = sqrt(accumSq);
  }

  // Loop over each singular vector
  for ( i=0; i<this.nRows; i++ ) {

    double wiInv = this.w[i];
    if ( NearZero(wiInv,1.0e-10) ) {
      continue;
    }

    double N = 0.0;           // N is the quasi-1-norm of the i-th column of this.U
    double alpha = 0.0;         // alpha is the dot product of this.dS and the i-th column of this.U

    const  *dSx = this.dS.GetPtr();
    const  *ux = this.U.GetColumnPtr(i);
    long j;
    for ( j=numEndEffectors; j>0; j-- ) {
      double tmp;
      alpha += (*ux)*(*(dSx++));
      tmp = Square( *(ux++) );
      alpha += (*ux)*(*(dSx++));
      tmp += Square(*(ux++));
      alpha += (*ux)*(*(dSx++));
      tmp += Square(*(ux++));
      N += sqrt(tmp);
    }

    // P is the quasi-1-norm of the response to angles changing according to the i-th column of this.V
    double P = 0.0;
    double *vx = this.V.GetColumnPtr(i);
    jnx = this.Jnorms.GetPtr();
    for ( j=this.nCols; j>0; j-- ) {
      double accum=0.0;
      for ( long k=numEndEffectors; k>0; k-- ) {
        accum += *(jnx++);
      }
      P += fabs((*(vx++)))*accum;
    }
  
    double lambda = 1.0;
    if ( N<P ) {
      lambda -= N/P;        // Scale back maximum permissable joint angle
    }
    lambda *= lambda;
    lambda *= this.DampingLambdaSDLS;

    // Calculate the this.dTheta from pure pseudoinverse considerations
    double scale = alpha*wiInv/(Square(wiInv)+Square(lambda));      // This times i-th column of this.V is the SDLS response
    MatrixRmn.AddArrayScale(this.nCols, this.V.GetColumnPtr(i), 1, this.dTheta.GetPtr(), 1, scale );
  }

  // Scale back to not exceed maximum angle changes
  double maxChange = this.dTheta.MaxAbs();
  if ( maxChange>MaxAngleSDLS ) {
    this.dTheta *= MaxAngleSDLS/maxChange;
  }
} */
}
