// @flow
import IKExample from './IKExample';
const THREE = require('three');

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
  75,
  window.innerWidth / window.innerHeight,
  0.1,
  1000
);

// const geometry = new THREE.BoxGeometry(1, 1, 1);
// const material = new THREE.MeshBasicMaterial({
//   color: (0x00ff00: number | string),
// });
// const cube = new THREE.Mesh(geometry, material);
// scene.add(cube);

camera.position.z = 5;

const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body && document.body.appendChild(renderer.domElement);

const ikExample = new IKExample(scene, 'IK_DLS');

function animate(timestamp: number) {
  requestAnimationFrame(animate);
  // todo: move this to timestep
  ikExample.stepSimulation(timestamp);
  ikExample.renderScene();
  renderer.render(scene, camera);
}
requestAnimationFrame(animate);
