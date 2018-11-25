// @flow

import {VectorR3} from './LinearR3';

type Purpose = 'JOINT' | 'EFFECTOR';

export const EFFECTOR = 'EFFECTOR';
export const JOINT = 'JOINT';

// class VectorR3;

export default class Node {
  // void PrintNode();
  // void InitNode();

  GetAttach() {
    return this.attach;
  }

  GetTheta() {
    return this.theta;
  }
  AddToTheta(delta: number) {
    this.theta += delta;
    return this.theta;
  }

  UpdateTheta(delta: number) {
    this.theta = delta;
    return this.theta;
  }

  GetS() {
    return this.s;
  }
  GetW() {
    return this.w;
  }

  GetMinTheta() {
    return this.minTheta;
  }
  GetMaxTheta() {
    return this.maxTheta;
  }
  GetRestAngle() {
    return this.restAngle;
  }
  SetTheta(newTheta: number) {
    this.theta = newTheta;
  }
  // void ComputeS(void);
  // void ComputeW(void);

  IsEffector() {
    return this.purpose === EFFECTOR;
  }
  IsJoint() {
    return this.purpose === JOINT;
  }
  GetEffectorNum() {
    return this.seqNumEffector;
  }
  GetJointNum() {
    return this.seqNumJoint;
  }

  IsFrozen() {
    return this.freezed;
  }
  Freeze() {
    this.freezed = true;
  }
  UnFreeze() {
    this.freezed = false;
  }

  //private:
  freezed: boolean; // Is this node frozen?
  seqNumJoint: number; // sequence number if this node is a joint
  seqNumEffector: number; // sequence number if this node is an effector
  size: number; // size
  purpose: Purpose; // joint / effector / both
  attach: VectorR3; // attachment point
  r: VectorR3; // relative position vector
  v: VectorR3; // rotation axis
  theta: number; // joint angle (radian)
  minTheta: number; // lower limit of joint angle
  maxTheta: number; // upper limit of joint angle
  restAngle: number; // rest position angle
  s: VectorR3; // GLobal Position
  w: VectorR3; // Global rotation axis
  left: ?Node; // left child
  right: ?Node; // right sibling
  realparent: ?Node; // pointer to real parent

  constructor(
    attach: VectorR3,
    v: VectorR3,
    size: number,
    purpose: Purpose,
    minTheta: number = -Math.PI,
    maxTheta: number = Math.PI,
    restAngle: number = 0
  ) {
    this.freezed = false;
    this.size = size;
    this.purpose = purpose;
    this.seqNumJoint = -1;
    this.seqNumEffector = -1;
    this.attach = new VectorR3(0.0, 0.0, 0.0).copy(attach); // Global attachment point when joints are at zero angle
    this.r = new VectorR3(0.0, 0.0, 0.0); // r will be updated when this node is inserted into tree
    this.v = new VectorR3(0.0, 0.0, 0.0).copy(v); // Rotation axis when joints at zero angles
    this.s = new VectorR3(0, 0, 0);
    this.w = new VectorR3(0, 0, 0);
    this.theta = 0.0;
    this.minTheta = minTheta;
    this.maxTheta = maxTheta;
    this.restAngle = restAngle;
    this.left = this.right = this.realparent = null;
  }

  // Compute the global position of a single node
  ComputeS() {
    let y: ?Node = this.realparent;
    let w: ?Node = this;
    this.s.copy(this.r); // Initialize to local (relative) position
    while (y) {
      this.s.Rotate(y.theta, y.v);
      y = y.realparent;
      w = w ? w.realparent : null;
      if (w != null) {
        this.s.add(w.r);
      }
    }
  }

  // Compute the global rotation axis of a single node
  ComputeW() {
    let y: ?Node = this.realparent;
    this.w.copy(this.v); // Initialize to local rotation axis
    while (y) {
      this.w.Rotate(y.theta, y.v);
      y = y.realparent;
    }
  }

  PrintNode() {
    console.log('Attach : (' + String(this.attach) + ')\n');
    console.log('r : (' + String(this.r) + ')\n');
    console.log('s : (' + String(this.s) + ')\n');
    console.log('w : (' + String(this.w) + ')\n');
    console.log(
      'realparent : ' +
        String(this.realparent ? this.realparent.seqNumJoint : 'no_parent') +
        '\n'
    );
  }

  InitNode() {
    this.theta = 0.0;
  }
}
