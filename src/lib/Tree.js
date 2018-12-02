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

import Node, {JOINT, EFFECTOR} from './Node';

function assert(pass: boolean) {
  if (!pass) {
    throw new Error('assertion failed');
  }
}

export default class Tree {
  // public:

  GetNumNode(): number {
    return this.nNode;
  }
  GetNumEffector(): number {
    return this.nEffector;
  }
  GetNumJoint(): number {
    return this.nJoint;
  }
  GetRoot() {
    return this.root;
  }
  GetParent(node: Node) {
    return node.realparent;
  }

  // private:
  root: ?Node = null;
  nNode: number = 0; // nNode = nEffector + nJoint
  nEffector: number = 0;
  nJoint: number = 0;

  // Initialize all nodes in the tree
  Init() {
    this.InitTree(this.root);
  }
  // Recursively initialize tree below the node
  InitTree(node: ?Node) {
    if (node != null) {
      node.InitNode();
      this.InitTree(node.left);
      this.InitTree(node.right);
    }
  }

  Compute() {
    this.ComputeTree(this.root);
  }

  ComputeTree(node: ?Node) {
    if (node != null) {
      node.ComputeS();
      node.ComputeW();
      this.ComputeTree(node.left);
      this.ComputeTree(node.right);
    }
  }

  InsertRoot(root: Node) {
    assert(this.nNode === 0);
    this.nNode++;
    this.root = root;
    root.r = root.attach.clone();
    assert(!(root.left || root.right));
    this.SetSeqNum(root);
  }

  InsertLeftChild(parent: Node, child: Node) {
    this.nNode++;
    parent.left = child;
    child.realparent = parent;
    child.r = child.attach.subtract(parent.attach);
    assert(!(child.left || child.right));
    this.SetSeqNum(child);
  }

  SetSeqNum(node: Node) {
    switch (node.purpose) {
      case JOINT:
        node.seqNumJoint = this.nJoint++;
        node.seqNumEffector = -1;
        break;
      case EFFECTOR:
        node.seqNumJoint = -1;
        node.seqNumEffector = this.nEffector++;
        break;
    }
  }

  GetSuccessor(fromNode: Node): ?Node {
    let node: ?Node = fromNode;
    if (node && node.left) {
      return node.left;
    }
    while (true) {
      if (node && node.right) {
        return node.right;
      }
      node = node ? node.realparent : null;
      if (!node) {
        return null; // Back to root, finished traversal
      }
    }
  }
}
