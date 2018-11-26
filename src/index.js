// @flow
import IKExample from './IKExample';
const THREE = require('three');
const OrbitControls = require('three-orbit-controls')(THREE);
const TransformControls = require('three-transform-controls')(THREE);

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
  45,
  window.innerWidth / window.innerHeight,
  1,
  10000
);

// const geometry = new THREE.BoxGeometry(1, 1, 1);
// const material = new THREE.MeshBasicMaterial({
//   color: (0x00ff00: number | string),
// });
// const cube = new THREE.Mesh(geometry, material);
// scene.add(cube);

camera.position.z = 5;

camera.position.set(10, 7, 10);
camera.lookAt(0, 0, 0);

const renderer = new THREE.WebGLRenderer();
renderer.setClearColor(0xeeeeee, 1);
renderer.setSize(window.innerWidth, window.innerHeight);
document.body && document.body.appendChild(renderer.domElement);
const debugTextEl = document.createElement('pre');
debugTextEl.id = 'debug';
document.body && document.body.appendChild(debugTextEl);

const transformControls = new TransformControls(camera, renderer.domElement);
const orbitControls = new OrbitControls(camera);
orbitControls.update();
scene.add(transformControls);

const debug = {log: 'hi'};
const ikExample = new IKExample('IK_SDLS', scene, transformControls, debug);

function animate(timestamp: number) {
  requestAnimationFrame(animate);
  // todo: move this to timestep
  ikExample.stepSimulation(timestamp);
  ikExample.renderScene();
  orbitControls.update();
  renderer.render(scene, camera);
  debugTextEl.textContent = debug.log;
  debug.log = '';
}
requestAnimationFrame(animate);

window.demo = {
  renderer,
  camera,
  scene,
  ikExample,
};
