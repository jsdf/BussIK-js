// @flow

import Node, {JOINT, EFFECTOR} from './bussik/Node';
import {
  VectorR3,
  VectorR3_UnitX,
  VectorR3_UnitY,
  VectorR3_UnitZ,
  VectorR3_Zero,
} from './bussik/LinearR3';

import Jacobian from './bussik/Jacobian';
import Tree from './bussik/Tree';

const THREE = require('three');

type IKMethod =
  | 'IK_JACOB_TRANS'
  | 'IK_PURE_PSEUDO'
  | 'IK_DLS'
  | 'IK_SDLS'
  | 'IK_DLS_SVD';

// declare class Tree {
//   Init(): void;
//   Compute(): void;
//   GetRoot(): Node;
//   InsertRoot(root: Node): void;
//   InsertLeftChild(parent: Node, child: Node): void;
// }

// declare class Node {
//   left: ?Node;
//   right: ?Node;
//   v: VectorR3;
//   r: VectorR3;
//   theta: number;

//   constructor(

//       new VectorR3(0.1, 0.0, 0.0875),
//       unitz,
//       0.08,
//       JOINT,
//       -1e30,
//       1e30,
//       radiansToDegrees(0)
//     )
// }

// declare class Jacobian {
//   constructor(tree: Tree): this;
//   Reset(): void;
//   SetJtargetActive(): void;
//   SetJendActive(): void;
//   ComputeJacobian(targetaa: Array<VectorR3>): void; // Set up Jacobian and deltaS vectors
//   CalcDeltaThetasTranspose(): void; // Jacobian transpose method
//   CalcDeltaThetasDLS(): void; // Damped least squares method
//   CalcDeltaThetasDLSwithSVD(): void;
//   CalcDeltaThetasPseudoinverse(): void; // Pure pseudoinverse method
//   CalcDeltaThetasSDLS(): void; // Selectively damped least squares method
//   ZeroDeltaThetas(): void;
//   UpdateThetas(): void; // Apply the change in the theta values
//   UpdatedSClampValue(targetaa: Array<VectorR3>): void;
// }

// declare class VectorR3 {
//   x: number;
//   y: number;
//   z: number;
//   constructor(x: number, y: number, z: number): this;
//   Set(x: number, y: number, z: number): void;
// }

function radiansToDegrees(angle: number) {
  return angle * (180 / Math.PI);
}

function nullthrows<T>(value: ?T): T {
  if (value == null) throw new Error('unexpected null');
  return value;
}

const MAX_NUM_NODE = 1000;
const MAX_NUM_THETA = 1000;
const MAX_NUM_EFFECT = 100;

let T = 0;
let targetaa: Array<VectorR3> = [];

// Make slowdown factor larger to make the simulation take larger, less frequent steps
// Make the constant factor in Tstep larger to make time pass more quickly
//const SlowdownFactor = 40;
const SlowdownFactor = 0; // Make higher to take larger steps less frequently
const SleepsPerStep = SlowdownFactor;
let SleepCounter = 0;
//const let Tstep = 0.0005*(let)SlowdownFactor;   // Time step

let AxesList: number; /* list to hold the axes    */
let AxesOn: boolean; /* ON or OFF        */

let Scale: number, Scale2: number; /* scaling factors      */

let JointLimitsOn: boolean;
let RestPositionOn: boolean;
let UseJacobianTargets1: boolean;

let numIteration: number = 1;
let error = 0.0;
let errorDLS = 0.0;
let errorSDLS = 0.0;
let sumError = 0.0;
let sumErrorDLS = 0.0;
let sumErrorSDLS = 0.0;

// #ifdef _DYNAMIC
// bool initMaxDist = true;
// extern let Excess[];
// extern let dsnorm[];
// #endif

function Reset(tree: Tree, ikJacobian: Jacobian) {
  AxesOn = false;

  Scale = 1.0;
  Scale2 = 0.0; /* because add 1. to it in Display()  */

  JointLimitsOn = true;
  RestPositionOn = false;
  UseJacobianTargets1 = false;

  tree.Init();
  tree.Compute();
  ikJacobian.Reset();
}

// Update target positions

function UpdateTargets(T2: number, treeY: Tree) {
  let T: number = T2 / 5;
  targetaa[0].Set(
    0.6 * Math.sin(0),
    0.6 * Math.cos(0),
    0.5 + 0.4 * Math.sin(3 * T)
  );
}

// Does a single update (on one kind of tree)
function DoUpdateStep(
  Tstep: number,
  treeY: Tree,
  jacob: Jacobian,
  ikMethod: IKMethod
) {
  if (SleepCounter == 0) {
    T += Tstep;
    UpdateTargets(T, treeY);
  }

  if (UseJacobianTargets1) {
    jacob.SetJtargetActive();
  } else {
    jacob.SetJendActive();
  }
  jacob.ComputeJacobian(targetaa); // Set up Jacobian and deltaS vectors

  // Calculate the change in theta values
  switch (ikMethod) {
    case 'IK_JACOB_TRANS':
      jacob.CalcDeltaThetasTranspose(); // Jacobian transpose method
      break;
    // case 'IK_DLS':
    //   jacob.CalcDeltaThetasDLS(); // Damped least squares method
    //   break;
    // case 'IK_DLS_SVD':
    //   jacob.CalcDeltaThetasDLSwithSVD();
    //   break;
    // case 'IK_PURE_PSEUDO':
    //   jacob.CalcDeltaThetasPseudoinverse(); // Pure pseudoinverse method
    //   break;
    // case 'IK_SDLS':
    //   jacob.CalcDeltaThetasSDLS(); // Selectively damped least squares method
    //   break;
    default:
      jacob.ZeroDeltaThetas();
      break;
  }

  if (SleepCounter == 0) {
    jacob.UpdateThetas(); // Apply the change in the theta values
    jacob.UpdatedSClampValue(targetaa);
    SleepCounter = SleepsPerStep;
  } else {
    SleepCounter--;
  }
}

function getLocalTransform(node: Node): THREE.Matrix4 {
  const transform = new THREE.Matrix4();
  let axis = new THREE.Vector3(node.v.x, node.v.y, node.v.z);
  let rotation = new THREE.Quaternion(0, 0, 0, 1);
  if (axis.length()) {
    transform.makeRotationAxis(axis, node.theta);
  } else {
    transform.makeRotationFromQuaternion(new THREE.Quaternion(0, 0, 0, 1));
  }

  transform.makeTranslation(node.r.x, node.r.y, node.r.z);
  return transform;
}

///quick demo showing the right-handed coordinate system and positive rotations around each axis
export default class IKExample {
  ikMethod: IKMethod;
  ikTree: Tree;
  ikNodes: Array<Node>;
  ikJacobian: Jacobian;

  targetInstance: THREE.Object3D;
  movingInstances: Array<THREE.Object3D>;
  scene: THREE.Scene;
  // public:

  constructor(scene: THREE.Scene, option: IKMethod) {
    this.scene = scene;
    this.ikMethod = option;

    ///create some graphics proxy for the tracking target
    ///the endeffector tries to track it using Inverse Kinematics
    this.targetInstance = this.makeTarget();

    this.BuildKukaIIWAShape();
    this.ikJacobian = new Jacobian(this.ikTree);
    Reset(this.ikTree, this.ikJacobian);
  }

  makeTarget(): THREE.Object3D {
    const geometry = new THREE.BoxGeometry(1, 1, 1);
    const material = new THREE.MeshBasicMaterial({
      color: (0x00ff00: number | string),
    });
    const cube = new THREE.Mesh(geometry, material);
    this.scene.add(cube);
    return cube;
  }

  MyDrawTree(node: ?Node, parentAccTransform: THREE.Matrix4): void {
    let lineColor = new VectorR3(0, 0, 0);
    let lineWidth = 2;
    if (node != null) {
      //  glPushMatrix();
      // const pos = new VectorR3(
      //   parentAccTransform.position.x,
      //   parentAccTransform.position.y,
      //   parentAccTransform.position.z
      // );

      // const color = new VectorR3(0, 1, 0);
      // const pointSize = 10;
      // this.app.this.renderer.drawPoint(pos, color, pointSize);

      // this.app.this.renderer.drawLine(
      //   pos,
      //   pos + 0.05 * parentAccTransform.getBasis().getColumn(0),
      //   new VectorR3(1, 0, 0),
      //   lineWidth
      // );
      // this.app.this.renderer.drawLine(
      //   pos,
      //   pos + 0.05 * parentAccTransform.getBasis().getColumn(1),
      //   new VectorR3(0, 1, 0),
      //   lineWidth
      // );
      // this.app.this.renderer.drawLine(
      //   pos,
      //   pos + 0.05 * parentAccTransform.getBasis().getColumn(2),
      //   new VectorR3(0, 0, 1),
      //   lineWidth
      // );

      // const axisLocal = new VectorR3(node.v.x, node.v.y, node.v.z);
      // const axisWorld = parentAccTransform.getBasis() * axisLocal;

      // this.app.this.renderer.drawLine(
      //   pos,
      //   pos + 0.1 * axisWorld,
      //   new VectorR3(0.2, 0.2, 0.7),
      //   5
      // );

      //node.DrawNode(node == root); // Recursively draw node and update ModelView matrix
      if (node.left) {
        const localTransform: THREE.Matrix4 = getLocalTransform(node.left);

        const leftAccTransform: THREE.Matrix4 = parentAccTransform.multiply(
          localTransform
        );
        // this.app.this.renderer.drawLine(
        //   parentAccTransform.getOrigin(),
        //   leftAccTransform.getOrigin(),
        //   lineColor,
        //   lineWidth
        // );
        this.MyDrawTree(node.left, leftAccTransform); // Draw tree of children recursively
      }
      //  glPopMatrix();
      if (node.right) {
        const localTransform = getLocalTransform(node.right);
        const rightAccTransform: THREE.Matrix4 = parentAccTransform.multiply(
          localTransform
        );
        // this.app.this.renderer.drawLine(
        //   parentAccTransform.getOrigin(),
        //   rightAccTransform.getOrigin(),
        //   lineColor,
        //   lineWidth
        // );
        this.MyDrawTree(node.right, rightAccTransform); // Draw right siblings recursively
      }
    }
  }
  stepSimulation(deltaTime: number): void {
    DoUpdateStep(deltaTime, this.ikTree, this.ikJacobian, this.ikMethod);
  }
  renderScene(): void {
    let localTransform: THREE.Matrix4 = getLocalTransform(
      nullthrows(this.ikTree.GetRoot())
    );
    this.MyDrawTree(this.ikTree.GetRoot(), localTransform);

    // let pos = new VectorR3(targetaa[0].x, targetaa[0].y, targetaa[0].z);
    // const orn: THREE.Quaternion = new THREE.Quaternion(0, 0, 0, 1);

    // this.app.this.renderer.writeSingleInstanceTransformToCPU(
    //   pos,
    //   orn,
    //   this.targetInstance
    // );
    // this.app.this.renderer.writeTransforms();
    // this.app.this.renderer.renderScene();
  }

  BuildKukaIIWAShape() {
    //const unitx:VectorR3  = VectorR3_UnitX;
    const unity: VectorR3 = VectorR3_UnitY;
    const unitz: VectorR3 = VectorR3_UnitZ;
    const unit1: VectorR3 = new VectorR3(
      Math.sqrt(14.0) / 8.0,
      1.0 / 8.0,
      7.0 / 8.0
    );
    const zero: VectorR3 = VectorR3_Zero;

    const minTheta = -4 * Math.PI;
    const maxTheta = 4 * Math.PI;

    // this.ikNodes.length = 8;//7DOF+additional endeffector

    this.ikNodes[0] = new Node(
      new VectorR3(0.1, 0.0, 0.0875),
      unitz,
      0.08,
      JOINT,
      -1e30,
      1e30,
      radiansToDegrees(0)
    );
    this.ikTree.InsertRoot(this.ikNodes[0]);

    this.ikNodes[1] = new Node(
      new VectorR3(0.1, -0.0, 0.29),
      unity,
      0.08,
      JOINT,
      -0.5,
      0.4,
      radiansToDegrees(0)
    );
    this.ikTree.InsertLeftChild(this.ikNodes[0], this.ikNodes[1]);

    this.ikNodes[2] = new Node(
      new VectorR3(0.1, -0.0, 0.4945),
      unitz,
      0.08,
      JOINT,
      minTheta,
      maxTheta,
      radiansToDegrees(0)
    );
    this.ikTree.InsertLeftChild(this.ikNodes[1], this.ikNodes[2]);

    this.ikNodes[3] = new Node(
      new VectorR3(0.1, 0.0, 0.71),
      unity.signFlip(),
      0.08,
      JOINT,
      minTheta,
      maxTheta,
      radiansToDegrees(0)
    );
    this.ikTree.InsertLeftChild(this.ikNodes[2], this.ikNodes[3]);

    this.ikNodes[4] = new Node(
      new VectorR3(0.1, 0.0, 0.8945),
      unitz,
      0.08,
      JOINT,
      minTheta,
      maxTheta,
      radiansToDegrees(0)
    );
    this.ikTree.InsertLeftChild(this.ikNodes[3], this.ikNodes[4]);

    this.ikNodes[5] = new Node(
      new VectorR3(0.1, 0.0, 1.11),
      unity,
      0.08,
      JOINT,
      minTheta,
      maxTheta,
      radiansToDegrees(0)
    );
    this.ikTree.InsertLeftChild(this.ikNodes[4], this.ikNodes[5]);

    this.ikNodes[6] = new Node(
      new VectorR3(0.1, 0.0, 1.191),
      unitz,
      0.08,
      JOINT,
      minTheta,
      maxTheta,
      radiansToDegrees(0)
    );
    this.ikTree.InsertLeftChild(this.ikNodes[5], this.ikNodes[6]);

    this.ikNodes[7] = new Node(
      new VectorR3(0.1, 0.0, 1.2),
      zero,
      0.08,
      EFFECTOR
    );
    this.ikTree.InsertLeftChild(this.ikNodes[6], this.ikNodes[7]);
  }
}
