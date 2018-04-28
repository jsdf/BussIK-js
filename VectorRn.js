// @flow

import type {VectorR3} from './LinearR3';

type Ptr = {
  postincrementDeref(): number,
};

export default class VectorRn extends Array<number> {
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
  SetTriple(i: number, vec: VectorR3) {
    this[i] = vec.x;
    this[i + 1] = vec.y;
    this[i + 2] = vec.z;
  }
  GetLength() {
    return this.length;
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
}
