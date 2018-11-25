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
// Linear Algebra Classes over R3
//
//
// A. Vector and Position classes
//
//    A.1. VectorR3: a real column vector of length 3.
//

// **************************************
// VectorR3 class                       *
// * * * * * * * * * * * * * * * * * * **

function Square(x) {
  return x * x;
}

export class VectorR3 {
  // public:
  // The x & y & z coordinates.
  x: number;
  y: number;
  z: number;
  constructor(x: number, y: number, z: number) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  toString() {
    return (
      '<' +
      this.x.toFixed(3) +
      ',' +
      this.y.toFixed(3) +
      ',' +
      this.z.toFixed(3) +
      '>'
    );
  }

  clone() {
    return new VectorR3(this.x, this.y, this.z);
  }

  copy(v: VectorR3) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;

    return this;
  }

  // VectorR3& operator+=
  add(v: VectorR3) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;
    return this;
  }

  // VectorR3& operator-=
  subtract(v: VectorR3) {
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;
    return this;
  }

  // VectorR3& operator*=
  multiplyScalar(m: number) {
    this.x *= m;
    this.y *= m;
    this.z *= m;
    return this;
  }

  // VectorR3 operator- () const { return ( VectorR3(-x, -y, -z) ); }
  signFlip() {
    return new VectorR3(-this.x, -this.y, -this.z);
  }
  // VectorR3& operator*= (const VectorR3& v); // Cross Product

  cross(v: VectorR3) {
    let tx = this.x;
    let ty = this.y;
    this.x = this.y * v.z - this.z * v.y;
    this.y = this.z * v.x - tx * v.z;
    this.z = tx * v.y - ty * v.x;
    return this;
  }

  Norm() {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
  }
  // double NormSq() const { return ( x*x + y*y + z*z ); }
  NormSq() {
    return this.x * this.x + this.y * this.y + this.z * this.z;
  }

  plus(u: VectorR3, v: VectorR3) {
    return new VectorR3(u.x + v.x, u.y + v.y, u.z + v.z);
  }

  // ******************************************************
  // * VectorR3 class - math library functions        *
  // * * * * * * * * * * * * * * * * * * * * * * * * * * **

  MaxAbs() {
    let m: number;
    m = this.x > 0.0 ? this.x : -this.x;
    if (this.y > m) m = this.y;
    else if (-this.y > m) m = -this.y;
    if (this.z > m) m = this.z;
    else if (-this.z > m) m = -this.z;
    return m;
  }

  Set(x: number, y: number, z: number) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  // *********************************************************************
  // Rotation routines                           *
  // *********************************************************************

  // s.Rotate(theta, u) rotates s and returns s
  //        rotated theta degrees around unit vector w.
  Rotate(theta: number, w: VectorR3) {
    let c: number = Math.cos(theta);
    let s: number = Math.sin(theta);
    let dotw: number = this.x * w.x + this.y * w.y + this.z * w.z;
    let v0x: number = dotw * w.x;
    let v0y: number = dotw * w.y; // v0 = provjection onto w
    let v0z: number = dotw * w.z;
    let v1x: number = this.x - v0x;
    let v1y: number = this.y - v0y; // v1 = projection onto plane normal to w
    let v1z: number = this.z - v0z;
    let v2x: number = w.y * v1z - w.z * v1y;
    let v2y: number = w.z * v1x - w.x * v1z; // v2 = w * v1 (cross product)
    let v2z: number = w.x * v1y - w.y * v1x;

    this.x = v0x + c * v1x + s * v2x;
    this.y = v0y + c * v1y + s * v2y;
    this.z = v0z + c * v1z + s * v2z;

    return this;
  }

  Dump(arraylike: Array<number>, offset: number) {
    arraylike[offset] = this.x;
    arraylike[offset + 1] = this.y;
    arraylike[offset + 2] = this.z;
  }
}

export const VectorR3_UnitVecIR3 = () => new VectorR3(1.0, 0.0, 0.0);
export const VectorR3_UnitVecJR3 = () => new VectorR3(0.0, 1.0, 0.0);
export const VectorR3_UnitVecKR3 = () => new VectorR3(0.0, 0.0, 1.0);

export const VectorR3_Zero = () => new VectorR3(0.0, 0.0, 0.0);
export const VectorR3_UnitX = () => new VectorR3(1.0, 0.0, 0.0);
export const VectorR3_UnitY = () => new VectorR3(0.0, 1.0, 0.0);
export const VectorR3_UnitZ = () => new VectorR3(0.0, 0.0, 1.0);
export const VectorR3_NegUnitX = () => new VectorR3(-1.0, 0.0, 0.0);
export const VectorR3_NegUnitY = () => new VectorR3(0.0, -1.0, 0.0);
export const VectorR3_NegUnitZ = () => new VectorR3(0.0, 0.0, -1.0);
