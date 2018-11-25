// @flow

import type {VectorR3} from './LinearR3';
import type MatrixRmn from './MatrixRmn';

type Ptr = {
  postincrementDeref(): number,
};

function assert(pass: boolean) {
  if (!pass) throw new Error('assertion failed');
}

let WorkVector: ?VectorRn = null;

export default class VectorRn extends Array<number> {
  // VectorRn::operator*=( double f )
  multiply(f: number) {
    for (var i = 0; i < this.length; i++) {
      this[i] *= f;
    }

    return this;
  }

  GetPtr(): Ptr {
    let cursor = 0;
    let array = this;
    return {
      postincrementDeref() {
        return array[cursor++];
      },
    };
  }
  SetLength(num: number) {
    if (this.length < num) {
      for (var i = this.length; i < num; i++) {
        this[i] = 0;
      }
    } else {
      this.length = num;
    }
  }
  // GetLength(): number;
  Fill(v: number) {
    this.fill(v);
  }
  SetZero() {
    this.fill(0);
  }
  SetTriple(i: number, u: VectorR3) {
    let j = 3 * i;
    assert(0 <= j && j + 2 < this.length);
    // u.Dump( x+j );
    this[j] = u.x;
    this[j + 1] = u.y;
    this[j + 2] = u.z;
  }
  GetLength() {
    return this.length;
  }

  LoadScaledMatrixRmn(D: MatrixRmn, d: number, scaleFactor: number) {
    let to /*this*/ = 0;
    for (let i = this.length; i > 0; i--) {
      this[to++] = D[d++] * scaleFactor;
    }
  }

  AddScaled(src: VectorRn, scaleFactor: number) {
    assert(src.length == this.length);
    let toPtr /*this*/ = 0;
    let fromPtr /*src*/ = 0;
    for (let i = this.length; i > 0; i--) {
      this[toPtr++] += src[fromPtr++] * scaleFactor;
    }
  }

  MaxAbs() {
    let result = 0.0;
    for (var i = 0; i < this.length; i++) {
      const abs = Math.abs(this[i]);
      if (abs > result) {
        result = abs;
      }
    }
    return result;
  }

  NormSq() {
    let res = 0.0;
    for (var i = 0; i < this.length; i++) {
      const val = this[i];
      res += val * val;
    }

    return res;
  }

  Norm() {
    return Math.sqrt(this.NormSq());
  }

  // STATIC METHODS

  static GetWorkVector(len: number) {
    if (WorkVector == null) {
      WorkVector = new VectorRn();
    }
    WorkVector.SetLength(len);
    return WorkVector;
  }
}
