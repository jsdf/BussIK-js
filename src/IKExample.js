// @flow

import Node, {JOINT, EFFECTOR} from './lib/Node';
import {
  VectorR3,
  VectorR3_UnitX,
  VectorR3_UnitY,
  VectorR3_UnitZ,
  VectorR3_Zero,
} from './lib/LinearR3';

import Jacobian from './lib/Jacobian';
import Tree from './lib/Tree';

const THREE = require('three');

type IKMethod =
  | 'IK_JACOB_TRANS'
  | 'IK_PURE_PSEUDO'
  | 'IK_DLS'
  | 'IK_SDLS'
  | 'IK_DLS_SVD';

function radiansToDegrees(angle: number) {
  return angle * (180 / Math.PI);
}

function nullthrows<T>(value: ?T): T {
  if (value == null) throw new Error('unexpected null');
  return value;
}

function debugPrintVector3(vec: VectorR3 | THREE.Vector3) {
  return `{x:${vec.x.toFixed(3)},y:${vec.y.toFixed(3)},z:${vec.z.toFixed(3)}}`;
}

const MAX_NUM_NODE = 1000;
const MAX_NUM_THETA = 1000;
const MAX_NUM_EFFECT = 100;

let T = 0;
let targetaa: Array<VectorR3> = [];
for (var i = 0; i < MAX_NUM_EFFECT; i++) {
  targetaa[i] = new VectorR3(0, 0, 0);
}

// Make slowdown factor larger to make the simulation take larger, less frequent steps
// Make the constant factor in Tstep larger to make time pass more quickly
//const SlowdownFactor = 40;
const SlowdownFactor = 2; // Make higher to take larger steps less frequently
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

const dummyRotation = new THREE.Quaternion();
const dummyScale = new THREE.Vector3();

function Reset(tree: Tree, ikJacobian: Jacobian) {
  AxesOn = false;

  Scale = 1.0;
  Scale2 = 0.0; /* because add 1. to it in Display()  */

  JointLimitsOn = true;
  RestPositionOn = false;
  UseJacobianTargets1 = true;

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
  // if (SleepCounter === 0) {
  //   T += Tstep;
  //   UpdateTargets(T, treeY);
  // }

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
    case 'IK_DLS':
      jacob.CalcDeltaThetasDLS(); // Damped least squares method
      break;
    // case 'IK_DLS_SVD':
    //   jacob.CalcDeltaThetasDLSwithSVD();
    //   break;
    // case 'IK_PURE_PSEUDO':
    //   jacob.CalcDeltaThetasPseudoinverse(); // Pure pseudoinverse method
    //   break;
    case 'IK_SDLS':
      jacob.CalcDeltaThetasSDLS(); // Selectively damped least squares method
      break;
    default:
      jacob.ZeroDeltaThetas();
      break;
  }

  if (SleepCounter === 0) {
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
  ikEndEffector: Node;
  ikJacobian: Jacobian;

  targetInstance: THREE.Object3D;
  nodesToLineVertices: Map<Node, THREE.Vector3> = new Map();
  lineGeometry: THREE.Geometry;
  debugPoints: Array<THREE.Object3D> = [];
  scene: THREE.Scene;
  transformControls: Object;
  debug: {log: string};

  // public:
  constructor(
    option: IKMethod,
    scene: THREE.Scene,
    transformControls: Object,
    debug: {log: string}
  ) {
    this.debug = debug;
    this.scene = scene;
    this.transformControls = transformControls;
    this.ikNodes = [];
    this.ikTree = new Tree();
    this.ikMethod = option;

    this.BuildKukaIIWAShape();
    this.ikJacobian = new Jacobian(this.ikTree);

    Reset(this.ikTree, this.ikJacobian);

    DoUpdateStep(1, this.ikTree, this.ikJacobian, this.ikMethod);
    this.createDebugPoints();
    this.updateDebugPoints();

    ///create some graphics proxy for the tracking target
    this.targetInstance = this.makeTarget(this.ikEndEffector.s);
  }

  createDebugPoints() {
    this.ikNodes.forEach(node => {
      const point = this.makeBox(0.1, Math.random() * 0xffffff);
      this.debugPoints.push(point);
      this.scene.add(point);
    });

    // build line
    const lineGeom = new THREE.Geometry();
    this.ikNodes.forEach(node => {
      const vert = new THREE.Vector3(0, 0, 0);
      lineGeom.vertices.push(vert);
    });
    const lineMaterial = new THREE.LineBasicMaterial({
      color: (0x0000ff: number | string),
    });

    const line = new THREE.Line(lineGeom, lineMaterial);
    this.scene.add(line);
    this.lineGeometry = lineGeom;
  }

  updateDebugPoints() {
    for (var i = 0; i < this.debugPoints.length; i++) {
      const pos = this.ikNodes[i].s;
      this.debugPoints[i].position.set(pos.x, pos.y, pos.z);
      this.lineGeometry.vertices[i].set(pos.x, pos.y, pos.z);
    }
    this.lineGeometry.verticesNeedUpdate = true;
  }

  makeBox(size: number, color: number | string, wireframe: boolean = false) {
    const geometry = new THREE.BoxGeometry(size, size, size);
    const material = new THREE.MeshBasicMaterial({color, wireframe});
    const cube = new THREE.Mesh(geometry, material);
    return cube;
  }

  makeTarget(pos: VectorR3): THREE.Object3D {
    const cube = this.makeBox(0.5, 0x00ff00, true);
    // cube.position.set(1, 1, 1);
    cube.position.set(pos.x, pos.y, pos.z);
    this.scene.add(cube);
    this.transformControls.attach(cube);
    return cube;
  }

  stepSimulation(deltaTime: number): void {
    targetaa[0].Set(
      this.targetInstance.position.x,
      this.targetInstance.position.y,
      this.targetInstance.position.z
    );
    // targetaa[0].Set(0.374, 0, 3);
    this.debug.log += `targetaa[0]${targetaa[0].toString()}\n`;
    DoUpdateStep(deltaTime, this.ikTree, this.ikJacobian, this.ikMethod);
  }

  renderScene(): void {
    this.debug.log += '\n';
    this.debug.log +=
      this.targetInstance &&
      `targetInstance.position: ${debugPrintVector3(
        this.targetInstance.position
      )}\n`;

    this.debug.log +=
      'ikNodes:\n' +
      this.ikNodes.map(node => node.s.toString()).join('\n') +
      '\n';

    this.updateDebugPoints();
    this.debug.log +=
      'debugPoints:\n' +
      this.debugPoints.map(p => debugPrintVector3(p.position)).join('\n') +
      '\n';
  }

  BuildKukaIIWAShape() {
    //const unitx:VectorR3  = VectorR3_UnitX();
    const unity: VectorR3 = VectorR3_UnitY();
    const unitz: VectorR3 = VectorR3_UnitZ();
    const unit1: VectorR3 = new VectorR3(
      Math.sqrt(14.0) / 8.0,
      1.0 / 8.0,
      7.0 / 8.0
    );
    const zero: VectorR3 = VectorR3_Zero();

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
    this.ikEndEffector = this.ikNodes[7];
  }
}
