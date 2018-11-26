// @flow
/*
*
* Mathematics Subpackage (VrMath)
*
*
* Author: Samuel R. Buss, sbuss@ucsd.edu.
* Web page: http://math.ucsd.edu/~sbuss/MathCG
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

//
// MatrixRmn:  Matrix over reals  (Variable dimensional vector)
//
//    Not very sophisticated yet.  Needs more functionality
//    To do: better handling of resizing.
//

import {VectorR3} from './LinearR3';

import VectorRn from './VectorRn';

const SVD_MAX_ITERATIONS = 100;

function assert(pass: boolean) {
  if (!pass) throw new Error('assertion failed');
}

function NearZero(x: number, tolerance: number) {
  return Math.abs(x) <= tolerance;
}

function Square(x: number) {
  return x * x;
}

interface Ptr {
  set(value: number): void;
  get(): number;
  increment(): Ptr;
  decrement(): Ptr;
  add(change: number): Ptr;
  subtract(change: number): Ptr;
  clone(): Ptr;
}

class ValuePtr implements Ptr {
  value: number;
  constructor(value) {
    this.value = value;
  }

  set(value) {
    this.value = value;
  }
  get() {
    return this.value;
  }

  increment() {
    throw new Error('unimplemented');
  }
  decrement() {
    throw new Error('unimplemented');
  }
  add(change: number) {
    throw new Error('unimplemented');
  }
  subtract(change: number) {
    throw new Error('unimplemented');
  }
  clone() {
    return new ValuePtr(this.value);
  }
}

class VectorRnPtr implements Ptr {
  _offset: number;
  _vec: VectorRn;
  constructor(vec: VectorRn, start: number = 0) {
    this._vec = vec;
    this._offset = start;
  }
  set(value) {
    this._vec[this._offset] = value;
  }
  get() {
    return this._vec[this._offset];
  }

  increment() {
    this._offset++;
    return this;
  }
  decrement() {
    this._offset--;
    return this;
  }
  add(change: number) {
    this._offset += change;
    return this;
  }
  subtract(change: number) {
    this._offset -= change;
    return this;
  }
  clone() {
    return new VectorRnPtr(this._vec, this._offset);
  }
}

let WorkMatrix = null;

function GetWorkMatrix(numRows: number, numCols: number) {
  if (WorkMatrix == null) {
    WorkMatrix = new MatrixRmn();
  }
  WorkMatrix.SetSize(numRows, numCols);
  return WorkMatrix;
}

// Helper routine to calculate dot product
function DotArray(
  length: number,
  A: MatrixRmn,
  ptrA: number,
  strideA: number,
  B: MatrixRmn,
  ptrB: number,
  strideB: number
) {
  let result = 0.0;
  for (; length > 0; length--) {
    result += A[ptrA] * B[ptrB];
    ptrA += strideA;
    ptrB += strideB;
  }
  return result;
}

// Helper routine: copies and scales an array (src and dest may be equal, or overlap)
function CopyArrayScale(
  length: number,
  fromM: MatrixRmn,
  fromPtr: number,
  fromStride: number,
  toM: MatrixRmn,
  toPtr: number,
  toStride: number,
  scale: number
) {
  for (; length > 0; length--) {
    toM[toPtr] = fromM[fromPtr] * scale;
    fromPtr += fromStride;
    toPtr += toStride;
  }
}

// Helper routine: adds a scaled array
//  fromArray = toArray*scale.
function AddArrayScale(
  length: number,
  fromM: MatrixRmn,
  fromPtr: number,
  fromStride: number,
  toM: MatrixRmn,
  toPtr: number,
  toStride: number,
  scale: number
) {
  for (; length > 0; length--) {
    toM[toPtr] += fromM[fromPtr] * scale;
    fromPtr += fromStride;
    toPtr += toStride;
  }
}

export default class MatrixRmn extends Array<number> {
  NumRows = 0;
  NumCols = 0;
  GetNumRows() {
    return this.NumRows;
  }
  GetNumColumns() {
    return this.NumCols;
  }

  GetColumnPtrOffset(j: number) {
    assert(0 <= j && j < this.NumCols);
    return j * this.NumRows;
  }

  SetSize(numRows: number, numCols: number) {
    assert(numRows > 0 && numCols > 0);
    this.length = numRows * numCols;
    this.NumRows = numRows;
    this.NumCols = numCols;
  }
  SetZero() {
    this.fill(0);
  }
  SetTriple(i: number, j: number, u: VectorR3): void {
    let ii = 3 * i;
    assert(0 <= i && ii + 2 < this.NumRows && 0 <= j && j < this.NumCols);
    // like this code (but working around the flow errors)
    // u.Dump(this, j * this.NumRows + ii);
    const start = j * this.NumRows + ii;
    this[start] = u.x;
    this[start + 1] = u.y;
    this[start + 2] = u.z;
  }

  SetIdentity() {
    assert(this.NumRows == this.NumCols);
    this.SetZero();
    this.SetDiagonalEntries(1.0);
  }

  // Fill the diagonal entries with the value d.  The rest of the matrix is unchanged.
  SetDiagonalEntries(d: number) {
    let diagLen = Math.min(this.NumRows, this.NumCols);
    let dPtr = 0; // this
    for (; diagLen > 0; diagLen--) {
      this[dPtr] = d;
      dPtr += this.NumRows + 1;
    }
  }

  // Multiply this matrix by column vector v.
  // Result is column vector "result"
  // void MatrixRmn::Multiply( const VectorRn& v, VectorRn& result ) const
  // {
  //   assert ( v.GetLength()==NumCols && result.GetLength()==NumRows );
  //   double* out = result.GetPtr();        // Points to entry in result vector
  //   const double* rowPtr = x;         // Points to beginning of next row in matrix
  //   for ( long j = NumRows; j>0; j-- ) {
  //     const double* in = v.GetPtr();
  //     const double* m = rowPtr++;
  //     *out = 0.0f;
  //     for ( long i = NumCols; i>0; i-- ) {
  //       *out += (*(in++)) * (*m);
  //       m += NumRows;
  //     }
  //     out++;
  //   }
  // }
  Multiply(v: VectorRn, result: VectorRn): void {
    assert(
      v.GetLength() === this.NumCols && result.GetLength() === this.NumRows
    );
    let resultPtr = 0; // Points to entry in result vector
    let rowPtr = 0; // Points to beginning of next row in matrix
    for (let j = this.NumRows; j > 0; j--) {
      let vPtr = 0; // ptr into v
      let m = rowPtr++; // ptr into this
      result[resultPtr] = 0.0;
      for (let i = this.NumCols; i > 0; i--) {
        result[resultPtr] += v[vPtr++] * this[m];
        m += this.NumRows;
      }
      resultPtr++;
    }
  }

  // Multiply transpose of this matrix by column vector v.
  //    Result is column vector "result"
  // Equivalent to mult by row vector on left
  // void MatrixRmn::MultiplyTranspose( const VectorRn& v, VectorRn& result ) const
  // {
  //   assert ( v.GetLength()==NumRows && result.GetLength()==NumCols );
  //   double* out = result.GetPtr();        // Points to entry in result vector
  //   const double* colPtr = x;         // Points to beginning of next column in matrix
  //   for ( long i=NumCols; i>0; i-- ) {
  //     const double* in=v.GetPtr();
  //     *out = 0.0f;
  //     for ( long j = NumRows; j>0; j-- ) {
  //       *out += (*(in++)) * (*(colPtr++));
  //     }
  //     out++;
  //   }
  // }
  MultiplyTranspose(v: VectorRn, result: VectorRn): void {
    assert(
      v.GetLength() === this.NumRows && result.GetLength() === this.NumCols
    );
    let resultPtr = 0; // Points to entry in result vector
    let colPtr = 0; // Points to beginning of next column in matrix
    for (let i = this.NumCols; i > 0; i--) {
      let vPtr = 0; // ptr into v
      result[resultPtr] = 0.0;
      for (let j = this.NumRows; j > 0; j--) {
        result[resultPtr] += v[vPtr++] * this[colPtr++];
      }
      resultPtr++;
    }
  }

  // The matrix A is loaded, in into "this" matrix, based at (0,0).
  //  The size of "this" matrix must be large enough to accomodate A.
  //  The rest of "this" matrix is left unchanged.  It is not filled with zeroes!

  LoadAsSubmatrix(A: MatrixRmn) {
    assert(A.NumRows <= this.NumRows && A.NumCols <= this.NumCols);
    let extraColStep = this.NumRows - A.NumRows;
    let toPtr = 0; // ptr into this
    let fromPtr = 0; // ptr into A
    for (let i = A.NumCols; i > 0; i--) {
      // Copy columns of A, one per time thru loop
      for (let j = A.NumRows; j > 0; j--) {
        // Copy all elements of this column of A
        this[toPtr++] = A[fromPtr++];
      }
      toPtr += extraColStep;
    }
  }

  // The matrix A is loaded, in transposed order into "this" matrix, based at (0,0).
  //  The size of "this" matrix must be large enough to accomodate A.
  //  The rest of "this" matrix is left unchanged.  It is not filled with zeroes!
  LoadAsSubmatrixTranspose(A: MatrixRmn) {
    assert(A.NumRows <= this.NumCols && A.NumCols <= this.NumRows);
    let rowPtr /*ptr into this*/ = 0;
    let fromPtr /*ptr into A*/ = 0;
    for (let i = A.NumCols; i > 0; i--) {
      // Copy columns of A, once per loop
      let toPtr /*ptr into this*/ = rowPtr;
      for (let j = A.NumRows; j > 0; j--) {
        // Loop copying values from the column of A
        this[toPtr] = A[fromPtr++];
        toPtr += this.NumRows;
      }
      rowPtr++;
    }
  }

  SetColumn(i: number, d: VectorRn) {
    assert(this.NumRows === d.GetLength());
    let toPtr = i * this.NumRows; // ptr into this
    let fromPtr = 0; // ptr into d
    for (i = this.NumRows; i > 0; i--) {
      this[toPtr++] = d[fromPtr++];
    }
  }

  // ConvertToRefNoFree
  // Converts the matrix (in place) to row echelon form
  // For us, row echelon form allows any non-zero values, not just 1's, in the
  //    position for a lead variable.
  // The "NoFree" version operates on the assumption that no free variable will be found.
  // Algorithm uses row operations and row pivoting (only).
  // Augmented matrix is correctly accomodated.  Only the first square part participates
  //    in the main work of row operations.
  ConvertToRefNoFree() {
    // Loop over all columns (variables)
    // Find row with most non-zero entry.
    // Swap to the highest active row
    // Subtract appropriately from all the lower rows (row op of type 3)
    let numIters = Math.min(this.NumRows, this.NumCols);
    let rowPtr1 = 0; // ptr into this
    let diagStep = this.NumRows + 1;
    let lenRowLeft = this.NumCols;
    for (; numIters > 1; numIters--) {
      // Find row with most non-zero entry.
      let rowPtr2 = rowPtr1; // ptr into this
      let maxAbs = Math.abs(this[rowPtr1]);
      let rowPivot = rowPtr1; // ptr into this
      let i;
      for (i = numIters - 1; i > 0; i--) {
        let newMax = this[++rowPivot];
        if (newMax > maxAbs) {
          maxAbs = this[rowPivot];
          rowPtr2 = rowPivot;
        } else if (-newMax > maxAbs) {
          maxAbs = -newMax;
          rowPtr2 = rowPivot;
        }
      }
      // Pivot step: Swap the row with highest entry to the current row
      if (rowPtr1 != rowPtr2) {
        let toPtr = rowPtr1; // ptr into this
        for (let i = lenRowLeft; i > 0; i--) {
          let temp = this[toPtr];
          this[toPtr] = this[rowPtr2];
          this[rowPtr2] = temp;
          toPtr += this.NumRows;
          rowPtr2 += this.NumRows;
        }
      }
      // Subtract this row appropriately from all the lower rows (row operation of type 3)
      rowPtr2 = rowPtr1;
      for (i = numIters - 1; i > 0; i--) {
        rowPtr2++;
        let toPtr = rowPtr2; // ptr into this
        let fromPtr = rowPtr1; // ptr into this
        assert(this[fromPtr] != 0.0);
        let alpha = this[toPtr] / this[fromPtr];
        this[toPtr] = 0.0;
        for (let j = lenRowLeft - 1; j > 0; j--) {
          toPtr += this.NumRows;
          fromPtr += this.NumRows;
          this[toPtr] -= this[fromPtr] * alpha;
        }
      }
      // Update for next iteration of loop
      rowPtr1 += diagStep;
      lenRowLeft--;
    }
  }

  // Adds d to each diagonal entry
  AddToDiagonal(d: number) {
    let diagLen = Math.min(this.NumRows, this.NumCols);
    let dPtr /*ptr into this*/ = 0;
    for (; diagLen > 0; diagLen--) {
      this[dPtr] += d;
      dPtr += this.NumRows + 1;
    }
    return this;
  }

  // Solves the equation   (*this)*xVec = b;
  // Uses row operations.  Assumes *this is square and invertible.
  // No error checking for divide by zero or instability (except with asserts)
  Solve(b: VectorRn, xVec: VectorRn) {
    assert(
      this.NumRows === this.NumCols &&
        this.NumCols === xVec.GetLength() &&
        this.NumRows === b.GetLength()
    );

    // Copy this matrix and b into an Augmented Matrix
    let AugMat = GetWorkMatrix(this.NumRows, this.NumCols + 1);
    AugMat.LoadAsSubmatrix(this);
    AugMat.SetColumn(this.NumRows, b);

    // Put into row echelon form with row operations
    AugMat.ConvertToRefNoFree();

    // Solve for x vector values using back substitution
    let xLast /* ptr into xVec */ = this.NumRows - 1; // Last entry in xVec
    let endRow /* ptr into AugMat */ = this.NumRows * this.NumCols - 1; // Last entry in the current row of the coefficient part of Augmented Matrix
    let bPtr /* ptr into AugMat */ = endRow + this.NumRows; // Last entry in augmented matrix (end of last column, in augmented part)
    for (let i = this.NumRows; i > 0; i--) {
      let accum = AugMat[bPtr--];
      // Next loop computes back substitution terms
      let rowPtr /* ptr into AugMat */ = endRow; // Points to entries of the current row for back substitution.
      let xPtr /* ptr into xVec */ = xLast; // Points to entries in the x vector (also for back substitution)
      for (let j = this.NumRows - i; j > 0; j--) {
        accum -= AugMat[rowPtr] * AugMat[xPtr--];
        rowPtr -= this.NumCols; // Previous entry in the row
      }
      assert(AugMat[rowPtr] != 0.0); // Are not supposed to be any free variables in this matrix
      xVec[xPtr] = accum / AugMat[rowPtr];
      endRow--;
    }
  }

  // Calculate the c=cosine and s=sine values for a Givens transformation.
  // The matrix M = ( (c, -s), (s, c) ) in row order transforms the
  //   column vector (a, b)^T to have y-coordinate zero.
  CalcGivensValues(a: number, b: number, c: Ptr, s: Ptr) {
    let denomInv = Math.sqrt(a * a + b * b);
    if (denomInv == 0.0) {
      c.set(1.0);
      s.set(0.0);
    } else {
      denomInv = 1.0 / denomInv;
      c.set(a * denomInv);
      s.set(-b * denomInv);
    }
  }

  // Applies Givens transform to columns i and i+1.
  // Equivalent to postmultiplying by the matrix
  //      ( c  -s )
  //    ( s   c )
  // with non-zero entries in rows i and i+1 and columns i and i+1
  PostApplyGivensAdj(c: number, s: number, idx: number) {
    assert(0 <= idx && idx < this.NumCols);
    let colA = idx * this.NumRows;
    let colB = colA + this.NumRows;
    for (let i = this.NumRows; i > 0; i--) {
      let temp = this[colA];
      this[colA] = this[colA] * c + this[colB] * s;
      this[colB] = this[colB] * c - temp * s;
      colA++;
      colB++;
    }
  }

  // Applies Givens transform to columns idx1 and idx2.
  // Equivalent to postmultiplying by the matrix
  //      ( c  -s )
  //    ( s   c )
  // with non-zero entries in rows idx1 and idx2 and columns idx1 and idx2
  PostApplyGivens(c: number, s: number, idx1: number, idx2: number) {
    assert(
      idx1 != idx2 &&
        0 <= idx1 &&
        idx1 < this.NumCols &&
        0 <= idx2 &&
        idx2 < this.NumCols
    );
    let colA = idx1 * this.NumRows;
    let colB = idx2 * this.NumRows;
    for (let i = this.NumRows; i > 0; i--) {
      let temp = this[colA];
      this[colA] = this[colA] * c + this[colB] * s;
      this[colB] = this[colB] * c - temp * s;
      colA++;
      colB++;
    }
  }

  // ********************************************************************************************
  // Singular value decomposition.
  // Return othogonal matrices U and V and diagonal matrix with diagonal w such that
  //     (this) = U * Diag(w) * V^T     (V^T is V-transpose.)
  // Diagonal entries have all non-zero entries before all zero entries, but are not
  //    necessarily sorted.  (Someday, I will write ComputedSortedSVD that handles
  //    sorting the eigenvalues by magnitude.)
  // ********************************************************************************************
  ComputeSVD(U: MatrixRmn, w: VectorRn, V: MatrixRmn) {
    assert(
      U.NumRows === this.NumRows &&
        V.NumCols === this.NumCols &&
        U.NumRows === U.NumCols &&
        V.NumRows === V.NumCols &&
        w.GetLength() === Math.min(this.NumRows, this.NumCols)
    );

    //  double temp=0.0;
    let superDiag = VectorRn.GetWorkVector(w.GetLength() - 1); // Some extra work space.  Will get passed around.

    // Choose larger of U, V to hold intermediate results
    // If U is larger than V, use U to store intermediate results
    // Otherwise use V.  In the latter case, we form the SVD of A transpose,
    //    (which is essentially identical to the SVD of A).
    let leftMatrix: MatrixRmn;
    let rightMatrix: MatrixRmn;
    if (this.NumRows >= this.NumCols) {
      U.LoadAsSubmatrix(this); // Copy A into U
      leftMatrix = U;
      rightMatrix = V;
    } else {
      V.LoadAsSubmatrixTranspose(this); // Copy A-transpose into V
      leftMatrix = V;
      rightMatrix = U;
    }

    // Do the actual work to calculate the SVD
    // Now matrix has at least as many rows as columns
    this.CalcBidiagonal(leftMatrix, rightMatrix, w, superDiag);
    this.ConvertBidiagToDiagonal(leftMatrix, rightMatrix, w, superDiag);
  }

  //TODO

  // ************************************************ CalcBidiagonal **************************
  // Helper routine for SVD computation
  // U is a matrix to be bidiagonalized.
  // On return, U and V are orthonormal and w holds the new diagonal
  //    elements and superDiag holds the super diagonal elements.

  CalcBidiagonal(U: MatrixRmn, V: MatrixRmn, w: VectorRn, superDiag: VectorRn) {
    assert(U.NumRows >= V.NumRows);

    // The diagonal and superdiagonal entries of the bidiagonalized
    //    version of the U matrix
    //    are stored in the vectors w and superDiag (temporarily).

    // Apply Householder transformations to U.
    // Householder transformations come in pairs.
    //   First, on the left, we map a portion of a column to zeros
    //   Second, on the right, we map a portion of a row to zeros
    const rowStep = U.NumCols;
    const diagStep = U.NumCols + 1;
    let diagPtr /*ptr into U*/ = 0;
    let wPtr /*ptr into w*/ = 0;
    let superDiagPtr /*ptr into superDiag*/ = 0;
    let colLengthLeft = U.NumRows;
    let rowLengthLeft = V.NumCols;
    while (true) {
      // Apply a Householder xform on left to zero part of a column
      this.SvdHouseholder(
        U,
        diagPtr,
        colLengthLeft,
        rowLengthLeft,
        1,
        rowStep,
        w,
        wPtr
      );

      if (rowLengthLeft == 2) {
        superDiag[superDiagPtr] = U[diagPtr + rowStep];
        break;
      }
      // Apply a Householder xform on the right to zero part of a row
      this.SvdHouseholder(
        U,
        diagPtr + rowStep,
        rowLengthLeft - 1,
        colLengthLeft,
        rowStep,
        1,
        superDiag,
        superDiagPtr
      );

      rowLengthLeft--;
      colLengthLeft--;
      diagPtr += diagStep;
      wPtr++;
      superDiagPtr++;
    }

    let extra = 0;
    diagPtr += diagStep;
    wPtr++;
    if (colLengthLeft > 2) {
      extra = 1;
      // Do one last Householder transformation when the matrix is not square
      colLengthLeft--;
      this.SvdHouseholder(U, diagPtr, colLengthLeft, 1, 1, 0, w, wPtr);
    } else {
      w[wPtr] = U[diagPtr];
    }

    // Form U and V from the Householder transformations
    V.ExpandHouseholders(V.NumCols - 2, 1, U, U.NumRows, U.NumRows, 1);
    U.ExpandHouseholders(V.NumCols - 1 + extra, 0, U, 0, 1, U.NumRows);

    // Done with bidiagonalization
    return;
  }

  // Helper routine for CalcBidiagonal
  // Performs a series of Householder transformations on a matrix
  // Stores results compactly into the matrix:   The Householder vector u (normalized)
  //   is stored into the first row/column being transformed.
  // The leading term of that row (= plus/minus its magnitude is returned
  //   separately into "retFirstEntry"
  SvdHouseholder(
    M1: MatrixRmn,
    basePt: number,
    colLength: number,
    numCols: number,
    colStride: number,
    rowStride: number,
    retV: VectorRn,
    retFirstEntry: number
  ) {
    // Calc norm of vector u
    let cPtr /*M1*/ = basePt;
    let norm = 0.0;
    let i;
    for (i = colLength; i > 0; i--) {
      norm += Square(M1[cPtr]);
      cPtr += colStride;
    }
    norm = Math.sqrt(norm); // Norm of vector to reflect to axis  e_1

    // Handle sign issues
    let imageVal; // Choose sign to maximize distance
    if (M1[basePt] < 0.0) {
      imageVal = norm;
      norm = 2.0 * norm * (norm - M1[basePt]);
    } else {
      imageVal = -norm;
      norm = 2.0 * norm * (norm + M1[basePt]);
    }
    norm = Math.sqrt(norm); // Norm is norm of reflection vector

    if (norm == 0.0) {
      // If the vector being transformed is equal to zero
      // Force to zero in case of roundoff errors
      cPtr = basePt;
      for (i = colLength; i > 0; i--) {
        M1[cPtr] = 0.0;
        cPtr += colStride;
      }
      retV[retFirstEntry] = 0.0;
      return;
    }

    retV[retFirstEntry] = imageVal;

    // Set up the normalized Householder vector
    M1[basePt] -= imageVal; // First component changes. Rest stay the same.
    // Normalize the vector
    norm = 1.0 / norm; // Now it is the inverse norm
    cPtr = basePt;
    for (i = colLength; i > 0; i--) {
      M1[cPtr] *= norm;
      cPtr += colStride;
    }

    // Transform the rest of the U matrix with the Householder transformation
    let rPtr /*M1*/ = basePt;
    for (let j = numCols - 1; j > 0; j--) {
      rPtr += rowStride;
      // Calc dot product with Householder transformation vector
      let dotP = DotArray(
        colLength,
        M1,
        basePt,
        colStride,
        M1,
        rPtr,
        colStride
      );
      // Transform with I - 2*dotP*(Householder vector)
      AddArrayScale(
        colLength,
        M1,
        basePt,
        colStride,
        M1,
        rPtr,
        colStride,
        -2.0 * dotP
      );
    }
  }

  // ********************************* ExpandHouseholders ********************************************
  // The matrix will be square.
  //   numXforms = number of Householder transformations to concatenate
  //    Each Householder transformation is represented by a unit vector
  //    Each successive Householder transformation starts one position later
  //      and has one more implied leading zero
  //   basePt = beginning of the first Householder transform
  //   colStride, rowStride: Householder xforms are stored in "columns"
  //   numZerosSkipped is the number of implicit zeros on the front each
  //      Householder transformation vector (only values supported are 0 and 1).
  ExpandHouseholders(
    numXforms: number,
    numZerosSkipped: number,
    M: MatrixRmn,
    basePt: number,
    colStride: number,
    rowStride: number
  ) {
    // Number of applications of the last Householder transform
    //     (That are not trivial!)
    let numToTransform = this.NumCols - numXforms + 1 - numZerosSkipped;
    assert(numToTransform > 0);

    if (numXforms == 0) {
      this.SetIdentity();
      return;
    }

    // Handle the first one separately as a special case,
    // "this" matrix will be treated to simulate being preloaded with the identity
    let hDiagStride = rowStride + colStride;
    let hBase /*M*/ = basePt + hDiagStride * (numXforms - 1); // Pointer to the last Householder vector
    let hDiagPtr /*M*/ = hBase + colStride * (numToTransform - 1); // Pointer to last entry in that vector
    let i;
    let diagPtr /*this*/ = this.NumCols * this.NumRows - 1; // Last entry in matrix (points to diagonal entry)
    let colPtr /*this*/ = diagPtr - (numToTransform - 1); // Pointer to column in matrix
    for (i = numToTransform; i > 0; i--) {
      CopyArrayScale(
        numToTransform,
        M,
        hBase,
        colStride,
        this,
        colPtr,
        1,
        -2.0 * M[hDiagPtr]
      );
      this[diagPtr] += 1.0; // Add back in 1 to the diagonal entry (since xforming the identity)
      diagPtr -= this.NumRows + 1; // Next diagonal entry in this matrix
      colPtr -= this.NumRows; // Next column in this matrix
      hDiagPtr -= colStride;
    }

    // Now handle the general case
    // A row of zeros must be in effect added to the top of each old column (in each loop)
    let colLastPtr /*this*/ = this.NumRows * this.NumCols - numToTransform - 1;
    for (i = numXforms - 1; i > 0; i--) {
      numToTransform++; // Number of non-trivial applications of this Householder transformation
      hBase -= hDiagStride; // Pointer to the beginning of the Householder transformation
      colPtr = colLastPtr;
      for (let j = numToTransform - 1; j > 0; j--) {
        // Get dot product
        let dotProd2N =
          -2.0 *
          DotArray(
            numToTransform - 1,
            M,
            hBase + colStride,
            colStride,
            this,
            colPtr + 1,
            1
          );
        this[colPtr] = dotProd2N * M[hBase]; // Adding onto zero at initial point
        AddArrayScale(
          numToTransform - 1,
          M,
          hBase + colStride,
          colStride,
          this,
          colPtr + 1,
          1,
          dotProd2N
        );
        colPtr -= this.NumRows;
      }
      // Do last one as a special case (may overwrite the Householder vector)
      CopyArrayScale(
        numToTransform,
        M,
        hBase,
        colStride,
        this,
        colPtr,
        1,
        -2.0 * M[hBase]
      );
      this[colPtr] += 1.0; // Add back one one as identity
      // Done with this Householder transformation
      colLastPtr--;
    }

    if (numZerosSkipped != 0) {
      assert(numZerosSkipped == 1);
      // Fill first row and column with identity (More generally: first numZerosSkipped many rows and columns)
      let d /*this*/ = 0;
      this[d] = 1;
      let d2 /*this*/ = d;
      for (i = this.NumRows - 1; i > 0; i--) {
        this[++d] = 0;
        this[(d2 += this.NumRows)] = 0;
      }
    }
  }

  // **************** ConvertBidiagToDiagonal ***********************************************
  // Do the iterative transformation from bidiagonal form to diagonal form using
  //    Givens transformation.  (Golub-Reinsch)
  // U and V are square.  Size of U less than or equal to that of U.
  ConvertBidiagToDiagonal(
    U: MatrixRmn,
    V: MatrixRmn,
    w: VectorRn,
    superDiag: VectorRn
  ) {
    // These two index into the last bidiagonal block  (last in the matrix, it will be
    //  first one handled.
    const state = {
      lastBidiagIdx: V.NumRows - 1,
      firstBidiagIdx: 0,
    };
    let eps = 1.0e-15 * Math.max(w.MaxAbs(), superDiag.MaxAbs());

    let iterations = 0; // infinite loops :(
    while (iterations < SVD_MAX_ITERATIONS) {
      iterations++;
      let workLeft = this.UpdateBidiagIndices(state, w, superDiag, eps);
      if (!workLeft) {
        break;
      }

      // Get ready for first Givens rotation
      // Push non-zero to M[2,1] with Givens transformation
      let wPtr = new VectorRnPtr(w, state.firstBidiagIdx);
      let sdPtr = new VectorRnPtr(superDiag, state.firstBidiagIdx);
      let extraOffDiag = new ValuePtr(0.0);
      if (wPtr.get() == 0.0) {
        this.ClearRowWithDiagonalZero(
          state.firstBidiagIdx,
          state.lastBidiagIdx,
          U,
          wPtr,
          sdPtr,
          eps
        );
        if (state.firstBidiagIdx > 0) {
          if (NearZero(sdPtr.decrement().get(), eps)) {
            sdPtr.set(0.0);
          } else {
            this.ClearColumnWithDiagonalZero(
              state.firstBidiagIdx,
              V,
              wPtr,
              sdPtr,
              eps
            );
          }
        }
        continue;
      }

      // Estimate an eigenvalue from bottom four entries of M
      // This gives a lambda value which will shift the Givens rotations
      // Last four entries of M^T * M are  ( ( A, B ), ( B, C ) ).
      let A;
      A =
        state.firstBidiagIdx < state.lastBidiagIdx - 1
          ? Square(superDiag[state.lastBidiagIdx - 2])
          : 0.0;
      let BSq = Square(w[state.lastBidiagIdx - 1]);
      A += BSq; // The "A" entry of M^T * M
      let C = Square(superDiag[state.lastBidiagIdx - 1]);
      BSq *= C; // The squared "B" entry
      C += Square(w[state.lastBidiagIdx]); // The "C" entry
      let lambda; // lambda will hold the estimated eigenvalue
      lambda = Math.sqrt(Square((A - C) * 0.5) + BSq); // Use the lambda value that is closest to C.
      if (A > C) {
        lambda = -lambda;
      }
      lambda += (A + C) * 0.5; // Now lambda equals the estimate for the last eigenvalue
      let t11 = Square(w[state.firstBidiagIdx]);
      let t12 = w[state.firstBidiagIdx] * superDiag[state.firstBidiagIdx];

      let c = new ValuePtr(0);
      let s = new ValuePtr(0);
      this.CalcGivensValues(t11 - lambda, t12, c, s);
      this.ApplyGivensCBTDCol(
        c.get(),
        s.get(),
        wPtr,
        sdPtr,
        extraOffDiag,
        wPtr.clone().add(1)
      );
      V.PostApplyGivensAdj(c.get(), -s.get(), state.firstBidiagIdx);
      let i;
      for (i = state.firstBidiagIdx; i < state.lastBidiagIdx - 1; i++) {
        // Push non-zero from M[i+1,i] to M[i,i+2]
        this.CalcGivensValues(wPtr.get(), extraOffDiag.get(), c, s);
        this.ApplyGivensCBTDRow(
          c.get(),
          s.get(),
          wPtr,
          sdPtr,
          extraOffDiag,
          extraOffDiag.get(),
          wPtr.clone().add(1),
          sdPtr.clone().add(1)
        );
        U.PostApplyGivensAdj(c.get(), -s.get(), i);
        // Push non-zero from M[i,i+2] to M[1+2,i+1]
        this.CalcGivensValues(sdPtr.get(), extraOffDiag.get(), c, s);
        this.ApplyGivensCBTDRow(
          c.get(),
          s.get(),
          sdPtr,
          wPtr.clone().add(1),
          extraOffDiag,
          extraOffDiag.get(),
          sdPtr.clone().add(1),
          wPtr.clone().add(2)
        );
        V.PostApplyGivensAdj(c.get(), -s.get(), i + 1);
        wPtr.increment();
        sdPtr.increment();
      }
      // Push non-zero value from M[i+1,i] to M[i,i+1] for i==lastBidiagIdx-1
      this.CalcGivensValues(wPtr.get(), extraOffDiag.get(), c, s);
      this.ApplyGivensCBTDCol(
        c.get(),
        s.get(),
        wPtr,
        extraOffDiag,
        sdPtr,
        wPtr.clone().add(1)
      );
      U.PostApplyGivensAdj(c.get(), -s.get(), i);

      // DEBUG
      // DebugCalcBidiagCheck( V, w, superDiag, U );
    }
  }

  // // This is called when there is a zero diagonal entry, with a non-zero superdiagonal entry on the same row.
  // // We use Givens rotations to "chase" the non-zero entry across the row; when it reaches the last
  // //  column, it is finally zeroed away.
  // // wPtr points to the zero entry on the diagonal.  sdPtr points to the non-zero superdiagonal entry on the same row.
  ClearRowWithDiagonalZero(
    firstBidiagIdx: number,
    lastBidiagIdx: number,
    U: MatrixRmn,
    wPtr: Ptr,
    sdPtr: Ptr,
    eps: number
  ) {
    let curSd = sdPtr.get(); // Value being chased across the row
    sdPtr.set(0.0);
    let i = firstBidiagIdx + 1;
    while (true) {
      // Rotate row i and row firstBidiagIdx (Givens rotation)
      let c = new ValuePtr(0);
      let s = new ValuePtr(0);
      this.CalcGivensValues(wPtr.increment().get(), curSd, c, s);
      U.PostApplyGivens(c.get(), -s.get(), i, firstBidiagIdx);
      wPtr.set(c.get() * wPtr.get() - s.get() * curSd);
      if (i == lastBidiagIdx) {
        break;
      }
      curSd = s.get() * sdPtr.increment().get(); // New value pops up one column over to the right
      sdPtr.set(c.get() * sdPtr.get());
      i++;
    }
  }

  // This is called when there is a zero diagonal entry, with a non-zero superdiagonal entry in the same column.
  // We use Givens rotations to "chase" the non-zero entry up the column; when it reaches the last
  //  column, it is finally zeroed away.
  // wPtr points to the zero entry on the diagonal.  sdPtr points to the non-zero superdiagonal entry in the same column.
  ClearColumnWithDiagonalZero(
    endIdx: number,
    V: MatrixRmn,
    wPtr: Ptr,
    sdPtr: Ptr,
    eps: number
  ) {
    let curSd = sdPtr.get(); // Value being chased up the column
    sdPtr.set(0.0);
    let i = endIdx - 1;
    while (true) {
      let c = new ValuePtr(0);
      let s = new ValuePtr(0);
      this.CalcGivensValues(wPtr.decrement().get(), curSd, c, s);
      V.PostApplyGivens(c.get(), -s.get(), i, endIdx);
      wPtr.set(c.get() * wPtr.get() - s.get() * curSd);
      if (i == 0) {
        break;
      }
      curSd = s.get() * sdPtr.decrement().get(); // New value pops up one row above
      if (NearZero(curSd, eps)) {
        break;
      }
      sdPtr.set(c.get() * sdPtr.get());
      i--;
    }
  }

  // // Matrix A is  ( ( a c ) ( b d ) ), i.e., given in column order.
  // // Mult's G[c,s]  times  A, replaces A.
  ApplyGivensCBTDCol(
    cosine: number,
    sine: number,
    a: Ptr,
    b: Ptr,
    c: Ptr,
    d: Ptr
  ) {
    let temp = a.get();
    a.set(cosine * a.get() - sine * b.get());
    b.set(sine * temp + cosine * b.get());
    temp = c.get();
    c.set(cosine * c.get() - sine * d.get());
    d.set(sine * temp + cosine * d.get());
  }

  // Now matrix A given in row order, A = ( ( a b c ) ( d e f ) ).
  // Return G[c,s] * A, replace A.  d becomes zero, no need to return.
  //  Also, it is certain the old *c value is taken to be zero!
  ApplyGivensCBTDRow(
    cosine: number,
    sine: number,
    a: Ptr,
    b: Ptr,
    c: Ptr,
    d: number,
    e: Ptr,
    f: Ptr
  ) {
    a.set(cosine * a.get() - sine * d);
    let temp = b.get();
    b.set(cosine * b.get() - sine * e.get());
    e.set(sine * temp + cosine * e.get());
    c.set(-sine * f.get());
    f.set(cosine * f.get());
  }

  // Helper routine for SVD conversion from bidiagonal to diagonal
  UpdateBidiagIndices(
    state: {firstBidiagIdx: number, lastBidiagIdx: number},
    w: VectorRn,
    superDiag: VectorRn,
    eps: number
  ) {
    let lastIdx = state.lastBidiagIdx;
    let sdPtr /*ptr into superDiag*/ = lastIdx - 1; // Entry above the last diagonal entry
    while (NearZero(superDiag[sdPtr], eps)) {
      superDiag[sdPtr--] = 0.0;
      lastIdx--;
      if (lastIdx == 0) {
        return false;
      }
    }
    state.lastBidiagIdx = lastIdx;
    let firstIdx = lastIdx - 1;
    let wPtr /*ptr into w*/ = firstIdx;
    while (firstIdx > 0) {
      if (NearZero(w[wPtr], eps)) {
        // If this diagonal entry (near) zero
        w[wPtr] = 0.0;
        break;
      }
      if (NearZero(superDiag[--sdPtr], eps)) {
        // If the entry above the diagonal entry is (near) zero
        superDiag[sdPtr] = 0.0;
        break;
      }
      wPtr--;
      firstIdx--;
    }
    state.firstBidiagIdx = firstIdx;
    return true;
  }

  // TODO
  // // ******************************************DEBUG STUFFF

  // bool MatrixRmn::DebugCheckSVD( const MatrixRmn& U, const VectorRn& w, const MatrixRmn& V ) const
  // {
  //   // Special SVD test code

  //   MatrixRmn IV( V.GetNumRows(), V.GetNumColumns() );
  //   IV.SetIdentity();
  //   MatrixRmn VTV( V.GetNumRows(), V.GetNumColumns() );
  //   MatrixRmn::TransposeMultiply( V, V, VTV );
  //   IV -= VTV;
  //   double error = IV.FrobeniusNorm();

  //   MatrixRmn IU( U.GetNumRows(), U.GetNumColumns() );
  //   IU.SetIdentity();
  //   MatrixRmn UTU( U.GetNumRows(), U.GetNumColumns() );
  //   MatrixRmn::TransposeMultiply( U, U, UTU );
  //   IU -= UTU;
  //   error += IU.FrobeniusNorm();

  //   MatrixRmn Diag( U.GetNumRows(), V.GetNumRows() );
  //   Diag.SetZero();
  //   Diag.SetDiagonalEntries( w );
  //   MatrixRmn B(U.GetNumRows(), V.GetNumRows() );
  //   MatrixRmn C(U.GetNumRows(), V.GetNumRows() );
  //   MatrixRmn::Multiply( U, Diag, B );
  //   MatrixRmn::MultiplyTranspose( B, V, C );
  //   C -= *this;
  //   error += C.FrobeniusNorm();

  //   bool ret = ( fabs(error)<=1.0e-13*w.MaxAbs() );
  //   assert ( ret );
  //   return ret;
  // }

  // bool MatrixRmn::DebugCheckInverse( const MatrixRmn& MInv ) const
  // {
  //     assert ( this->NumRows==this->NumCols );
  //     assert ( MInv.NumRows==MInv.NumCols );
  //     MatrixRmn I(this->NumRows, this->NumCols);
  //     I.SetIdentity();
  //     MatrixRmn MMInv(this->NumRows, this->NumCols);
  //     Multiply(*this, MInv, MMInv);
  //     I -= MMInv;
  //     double error = I.FrobeniusNorm();
  //     bool ret = ( fabs(error)<=1.0e-13 );
  //     assert ( ret );
  //     return ret;
  // }

  // bool MatrixRmn::DebugCalcBidiagCheck( const MatrixRmn& U, const VectorRn& w, const VectorRn& superDiag, const MatrixRmn& V ) const
  // {
  //   // Special SVD test code

  //   MatrixRmn IV( V.GetNumRows(), V.GetNumColumns() );
  //   IV.SetIdentity();
  //   MatrixRmn VTV( V.GetNumRows(), V.GetNumColumns() );
  //   MatrixRmn::TransposeMultiply( V, V, VTV );
  //   IV -= VTV;
  //   double error = IV.FrobeniusNorm();

  //   MatrixRmn IU( U.GetNumRows(), U.GetNumColumns() );
  //   IU.SetIdentity();
  //   MatrixRmn UTU( U.GetNumRows(), U.GetNumColumns() );
  //   MatrixRmn::TransposeMultiply( U, U, UTU );
  //   IU -= UTU;
  //   error += IU.FrobeniusNorm();

  //   MatrixRmn DiagAndSuper( U.GetNumRows(), V.GetNumRows() );
  //   DiagAndSuper.SetZero();
  //   DiagAndSuper.SetDiagonalEntries( w );
  //   if ( this->GetNumRows()>=this->GetNumColumns() ) {
  //     DiagAndSuper.SetSequence( superDiag, 0, 1, 1, 1 );
  //   }
  //   else {
  //     DiagAndSuper.SetSequence( superDiag, 1, 0, 1, 1 );
  //   }
  //   MatrixRmn B(U.GetNumRows(), V.GetNumRows() );
  //   MatrixRmn C(U.GetNumRows(), V.GetNumRows() );
  //   MatrixRmn::Multiply( U, DiagAndSuper, B );
  //   MatrixRmn::MultiplyTranspose( B, V, C );
  //   C -= *this;
  //   error += C.FrobeniusNorm();

  //   bool ret = ( fabs(error)<1.0e-13*Max(w.MaxAbs(),superDiag.MaxAbs()) );
  //   assert ( ret );
  //   return ret;
  // }

  // STATIC METHODS

  // Multiply two MatrixRmn's.  Transpose the second matrix before multiplying
  static MultiplyTranspose(A: MatrixRmn, B: MatrixRmn, dst: MatrixRmn) {
    assert(
      A.NumCols == B.NumCols &&
        A.NumRows == dst.NumRows &&
        B.NumRows == dst.NumCols
    );
    let length = A.NumCols;

    let bPtr /*B*/ = 0; // Points to beginning of row in B
    let dPtr /*dst*/ = 0;
    for (let i = dst.NumCols; i > 0; i--) {
      let aPtr /*A*/ = 0; // Points to beginning of row in A
      for (let j = dst.NumRows; j > 0; j--) {
        dst[dPtr] = DotArray(length, A, aPtr, A.NumRows, B, bPtr, B.NumRows);
        dPtr++;
        aPtr++;
      }
      bPtr++;
    }

    return dst;
  }
}
