'use strict';

let gl;                         // The webgl context.
let surface;                    // A surface model
let shProgram;                  // A shader program
let spaceball;                  // A SimpleRotator object that lets the user rotate the view by mouse.
let video;
let track;
let texture1, texture2;
let virtualCam, webCam;
function deg2rad(angle) {
  return angle * Math.PI / 180;
}


// Constructor
function Model(name) {
  this.name = name;
  this.iVertexBuffer = gl.createBuffer();
  this.iTextureBuffer = gl.createBuffer();
  this.countT = 0;
  this.count = 0;

  this.BufferData = function(vertices) {

    gl.bindBuffer(gl.ARRAY_BUFFER, this.iVertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STREAM_DRAW);

    this.count = vertices.length / 3;
  }

  this.TextureBufferData = function(points) {

    gl.bindBuffer(gl.ARRAY_BUFFER, this.iTextureBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(points), gl.STREAM_DRAW);

    this.countT = points.length / 2;
  }

  this.Draw = function() {

    gl.bindBuffer(gl.ARRAY_BUFFER, this.iVertexBuffer);
    gl.vertexAttribPointer(shProgram.iAttribVertex, 3, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(shProgram.iAttribVertex);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.iTextureBuffer);
    gl.vertexAttribPointer(shProgram.iAttribTexture, 2, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(shProgram.iAttribTexture);

    gl.drawArrays(gl.TRIANGLE_STRIP, 0, this.count);
  }

  this.DrawPoint = function() {
    gl.bindBuffer(gl.ARRAY_BUFFER, this.iVertexBuffer);
    gl.vertexAttribPointer(shProgram.iAttribVertex, 3, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(shProgram.iAttribVertex);
    gl.drawArrays(gl.LINE_STRIP, 0, this.count);
  }
}


// Constructor
function ShaderProgram(name, program) {

  this.name = name;
  this.prog = program;

  // Location of the attribute variable in the shader program.
  this.iAttribVertex = -1;
  this.iAttribTexture = -1;
  // Location of the uniform matrix representing the combined transformation.
  this.iModelViewProjectionMatrix = -1;

  this.iTMU = -1;
  this.iScale = 1.0;

  this.Use = function() {
    gl.useProgram(this.prog);
  }
}


/* Draws a colored cube, along with a set of coordinate axes.
 * (Note that the use of the above drawPrimitive function is not an efficient
 * way to draw with WebGL.  Here, the geometry is so simple that it doesn't matter.)
 */
function draw() {
  let spanValues = document.getElementsByClassName("spanValue");
  let eyeSeparation = 70.0;
  eyeSeparation = document.getElementById("ES").value;
  spanValues[0].innerHTML = eyeSeparation;
  virtualCam.mEyeSeparation = eyeSeparation;
  let ratio = 1.0;
  let convergence = 2000.0;
  convergence = document.getElementById("C").value;
  spanValues[1].innerHTML = convergence;
  virtualCam.mConvergence = convergence

  let nearClippingDistance = 5.0;
  nearClippingDistance = document.getElementById("NCD").value;
  spanValues[2].innerHTML = nearClippingDistance;
  console.log(nearClippingDistance)
  virtualCam.mNearClippingDistance = parseFloat(nearClippingDistance);
  let fov = 0.8;
  fov = document.getElementById("FOV").value;
  spanValues[3].innerHTML = fov;
  virtualCam.mFOV = fov;

  gl.clearColor(0, 0, 0, 1);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

  /* Set the values of the projection transformation */
  let para = 3
  let projection = m4.orthographic(-para, para, -para, para, 0, para * 4);
  let tSS = vectorFromVec(vec)

  /* Get the view matrix from the SimpleRotator object.*/
  let modelView = spaceball.getViewMatrix();
  let initMatrix = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]

  let rotateToPointZero = m4.axisRotation([0.707, 0.707, 0], 0.);
  let translateToPointLeft = m4.translation(-3, -3, -10);
  let translateToPointZero = m4.translation(0, 0, -11.5);
  let matAccum0 = m4.multiply(rotateToPointZero, initMatrix);
  /*if (matrix != null) {
     matAccum0 = m4.multiply(rotateToPointZero, matrix);
   }*/
  let matAccum1 = m4.multiply(translateToPointLeft, matAccum0);
  let matAccum2 = m4.multiply(translateToPointZero, matAccum0);
  let matAccum3 = m4.multiply(m4.multiply(m4.translation(tSS[0], tSS[1], tSS[2]), translateToPointZero), matAccum0);

  /* Multiply the projection matrix times the modelview matrix to give the
     combined transformation matrix, and send that to the shader program. */
  gl.uniform1i(shProgram.iTMU, 0);
  gl.enable(gl.TEXTURE_2D);
  let modelViewProjection = m4.multiply(projection, matAccum1);
  gl.uniformMatrix4fv(shProgram.iModelViewProjectionMatrix, false, modelViewProjection);
  
  gl.bindTexture(gl.TEXTURE_2D, texture1);
  virtualCam.ApplyLeftFrustum();
  gl.uniformMatrix4fv(shProgram.iModelViewProjectionMatrix, false, m4.multiply(modelView, m4.multiply(virtualCam.mProjectionMatrix, matAccum2)));
  gl.colorMask(false, true, true, false);
  surface.Draw();
  gl.clear(gl.DEPTH_BUFFER_BIT);
  virtualCam.ApplyRightFrustum();
  let finalMat = m4.multiply(modelView, m4.multiply(virtualCam.mProjectionMatrix, matAccum2))
  let finalMatSphere = m4.multiply(modelView, m4.multiply(virtualCam.mProjectionMatrix, matAccum3))
  gl.uniformMatrix4fv(shProgram.iModelViewProjectionMatrix, false, finalMat);
  gl.colorMask(true, false, false, false);
  surface.Draw();
  gl.colorMask(true, true, true, true);

  gl.uniformMatrix4fv(shProgram.iModelViewProjectionMatrix, false, finalMatSphere);
  sph.Draw();
  if (panner) {
    panner.setPosition(tSS[0] * 1.5, tSS[1] * 1.5, tSS[2] * 1.5);
  }
  
}

function animate() {
  draw()
  window.requestAnimationFrame(animate)
}

function dot(a, b) {
  let c = [(a[1] * b[2] - a[2] * b[1]), (a[0] * b[2] - b[0] * a[2]), (a[0] * b[1] - a[1] * b[0])]
  return c
}
function normalize(a) {
  let d = Math.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)
  let n = [a[0] / d, a[1] / d, a[2] / d]
  return n;
}

function CreateSurfaceData() {
  let vertexList = [];
  let i = 0;
  let j = 0;
  let a = -0.2633257397764612;
  let c = 3.2;
  let b = 1.6099263856487789;
  while (i < b) {
    while (j < Math.PI * 2) {
      let v1 = conjugation(i, j, a, c)
      let v2 = conjugation(i + 0.2, j, a, c)
      let v3 = conjugation(i, j + 0.2, a, c)
      let v4 = conjugation(i + 0.2, j + 0.2, a, c)
      vertexList.push(v1.x, v1.y, v1.z);
      vertexList.push(v2.x, v2.y, v2.z);
      vertexList.push(v3.x, v3.y, v3.z);
      vertexList.push(v3.x, v3.y, v3.z);
      vertexList.push(v4.x, v4.y, v4.z);
      vertexList.push(v2.x, v2.y, v2.z);
      j += 0.2
    }
    j = 0;
    i += 0.2
  }
  return vertexList;
}
function CreateTextureData() {
  let texCoordList = [];
  let i = 0;
  let j = 0;
  let a = -0.2633257397764612;
  let c = 3.2;
  let b = 1.6099263856487789;
  while (i < b) {
    while (j < Math.PI * 2) {
      let u = map(i, 0, b, 0, 1);
      let v = map(j, 0, Math.PI * 2, 0, 1);
      texCoordList.push(u, v);
      u = map(i + 0.2, 0, b, 0, 1);
      texCoordList.push(u, v);
      u = map(i, 0, b, 0, 1);
      v = map(j + 0.2, 0, Math.PI * 2, 0, 1);
      texCoordList.push(u, v);
      texCoordList.push(u, v);
      u = map(i + 0.2, 0, b, 0, 1);
      v = map(j + 0.2, 0, Math.PI * 2, 0, 1);
      texCoordList.push(u, v);
      u = map(i + 0.2, 0, b, 0, 1);
      v = map(j, 0, Math.PI * 2, 0, 1);
      texCoordList.push(u, v);

      j += 0.2;
    }
    j = 0
    i += 0.2;
  }
  return texCoordList;
}

function conjugation(z, b, a, c) {

  let r = a * (1 - Math.cos(Math.PI * 2 * z / c)) + 1;
  let t = 0.8 * Math.PI;
  let x = 0.8 * r * Math.cos(b);
  let y = 0.8 * r * Math.sin(b)
  let z1 = 0.8 * z
  return { x: x, y: y, z: z1 }
}

function map(val, f1, t1, f2, t2) {
  let m;
  m = (val - f1) * (t2 - f2) / (t1 - f1) + f2
  return Math.min(Math.max(m, f2), t2);
}

function vec3Cross(a, b) {
  let x = a.y * b.z - b.y * a.z;
  let y = a.z * b.x - b.z * a.x;
  let z = a.x * b.y - b.x * a.y;
  return { x: x, y: y, z: z }
}

function vec3Normalize(a) {
  var mag = Math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  a[0] /= mag; a[1] /= mag; a[2] /= mag;
}

/* Initialize the WebGL context. Called from init() */
let sph;
function initGL() {
  let prog = createProgram(gl, vertexShaderSource, fragmentShaderSource);

  shProgram = new ShaderProgram('Basic', prog);
  shProgram.Use();

  shProgram.iAttribVertex = gl.getAttribLocation(prog, "vertex");
  shProgram.iAttribTexture = gl.getAttribLocation(prog, "texture");
  shProgram.iModelViewProjectionMatrix = gl.getUniformLocation(prog, "ModelViewProjectionMatrix");
  shProgram.iTMU = gl.getUniformLocation(prog, 'tmu');
  shProgram.iScale = gl.getUniformLocation(prog, 'scl');

  virtualCam = new StereoCamera(2000, 70.0, 1, 0.8, 5, 15);

  surface = new Model('Surface');
  surface.BufferData(CreateSurfaceData());
  LoadTexture();
  surface.TextureBufferData(CreateTextureData());
  gl.enable(gl.DEPTH_TEST);

  webCam = new Model("Web Camera")

  let siz = 6;
  webCam.BufferData([0, 0, 0, siz, 0, 0, siz, siz, 0, siz, siz, 0, 0, siz, 0, 0, 0, 0]);
  webCam.TextureBufferData([1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1]);
  sph = new Model("Sound Source");
  sph.BufferData(createSphereSurface(0.1));
  sph.TextureBufferData(createSphereSurface(0.1));
}
function createSphereSurface(r) {
  let vertexList = [];
  let lon = -Math.PI;
  let lat = -Math.PI * 0.5;
  const STEP = 0.1;
  while (lon < Math.PI) {
    while (lat < Math.PI * 0.5) {
      let v1 = sphere(r, lon, lat);
      let v2 = sphere(r, lon + STEP, lat);
      let v3 = sphere(r, lon, lat + STEP);
      let v4 = sphere(r, lon + STEP, lat + STEP);
      vertexList.push(v1.x, v1.y, v1.z);
      vertexList.push(v2.x, v2.y, v2.z);
      vertexList.push(v3.x, v3.y, v3.z);
      vertexList.push(v3.x, v3.y, v3.z);
      vertexList.push(v4.x, v4.y, v4.z);
      vertexList.push(v2.x, v2.y, v2.z);
      lat += STEP;
    }
    lat = -Math.PI * 0.5
    lon += STEP;
  }
  return vertexList;
}

function sphere(r, u, v) {
  let x = r * Math.sin(u) * Math.cos(v);
  let y = r * Math.sin(u) * Math.sin(v);
  let z = r * Math.cos(u);
  return { x: x, y: y, z: z };
}


/* Creates a program for use in the WebGL context gl, and returns the
 * identifier for that program.  If an error occurs while compiling or
 * linking the program, an exception of type Error is thrown.  The error
 * string contains the compilation or linking error.  If no error occurs,
 * the program identifier is the return value of the function.
 * The second and third parameters are strings that contain the
 * source code for the vertex shader and for the fragment shader.
 */
function createProgram(gl, vShader, fShader) {
  let vsh = gl.createShader(gl.VERTEX_SHADER);
  gl.shaderSource(vsh, vShader);
  gl.compileShader(vsh);
  if (!gl.getShaderParameter(vsh, gl.COMPILE_STATUS)) {
    throw new Error("Error in vertex shader:  " + gl.getShaderInfoLog(vsh));
  }
  let fsh = gl.createShader(gl.FRAGMENT_SHADER);
  gl.shaderSource(fsh, fShader);
  gl.compileShader(fsh);
  if (!gl.getShaderParameter(fsh, gl.COMPILE_STATUS)) {
    throw new Error("Error in fragment shader:  " + gl.getShaderInfoLog(fsh));
  }
  let prog = gl.createProgram();
  gl.attachShader(prog, vsh);
  gl.attachShader(prog, fsh);
  gl.linkProgram(prog);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
    throw new Error("Link error in program:  " + gl.getProgramInfoLog(prog));
  }
  return prog;
}


/**
 * initialization function that will be called when the page has loaded
 */
function init() {
  let canvas;
  try {
    let resolution = Math.min(window.innerHeight, window.innerWidth);
    canvas = document.querySelector('canvas');
    gl = canvas.getContext("webgl");
    playSomeMusic()
    gl.viewport(0, 0, 600, 600);
    if (!gl) {
      throw "Browser does not support WebGL";
    }
  }
  catch (e) {
    document.getElementById("canvas-holder").innerHTML = e
    return;
  }
  try {
    initGL();  // initialize the WebGL graphics context
  }
  catch (e) {
    document.getElementById("canvas-holder").innerHTML =
      "<p>Sorry, could not initialize the WebGL graphics context: " + e + "</p>";
    return;
  }

  spaceball = new TrackballRotator(canvas, draw, 0);

  animate()
}

function mat4Transpose(a, transposed) {
  var t = 0;
  for (var i = 0; i < 4; ++i) {
    for (var j = 0; j < 4; ++j) {
      transposed[t++] = a[j * 4 + i];
    }
  }
}

function mat4Invert(m, inverse) {
  var inv = new Float32Array(16);
  inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] +
    m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
  inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] -
    m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
  inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] +
    m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
  inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] -
    m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
  inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] -
    m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
  inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] +
    m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
  inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] -
    m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
  inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] +
    m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
  inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] +
    m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
  inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] -
    m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
  inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] +
    m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
  inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] -
    m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
  inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] -
    m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
  inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] +
    m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
  inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] -
    m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
  inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] +
    m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

  var det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
  if (det == 0) return false;
  det = 1.0 / det;
  for (var i = 0; i < 16; i++) inverse[i] = inv[i] * det;
  return true;
}

function LoadTexture() {
  texture1 = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, texture1);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);

  const image = new Image();
  image.crossOrigin = 'anonymus';
  image.src =  "https://raw.githubusercontent.com/anton-taranec/GW/main/texture.jpg";
  image.onload = () => {
    gl.bindTexture(gl.TEXTURE_2D, texture1);
    gl.texImage2D(
      gl.TEXTURE_2D,
      0,
      gl.RGBA,
      gl.RGBA,
      gl.UNSIGNED_BYTE,
      image
    );
    draw()
  }
}

function getWebcam() {
  navigator.getUserMedia({ video: true, audio: false }, function(stream) {
    video.srcObject = stream;
    track = stream.getTracks()[0];
  }, function(e) {
    console.error('Rejected!', e);
  });
}

function CreateWebCamTexture() {
  texture2 = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, texture2);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  return texture2;
}

function StereoCamera(
  Convergence,
  EyeSeparation,
  AspectRatio,
  FOV,
  NearClippingDistance,
  FarClippingDistance
) {
  this.mConvergence = Convergence;
  this.mEyeSeparation = EyeSeparation;
  this.mAspectRatio = AspectRatio;
  this.mFOV = FOV;
  this.mNearClippingDistance = NearClippingDistance;
  this.mFarClippingDistance = FarClippingDistance;

  this.mProjectionMatrix = null;
  this.mModelViewMatrix = null;

  this.ApplyLeftFrustum = function() {
    let top, bottom, left, right;
    top = this.mNearClippingDistance * Math.tan(this.mFOV / 2);
    bottom = -top;

    let a = this.mAspectRatio * Math.tan(this.mFOV / 2) * this.mConvergence;
    let b = a - this.mEyeSeparation / 2;
    let c = a + this.mEyeSeparation / 2;

    left = (-b * this.mNearClippingDistance) / this.mConvergence;
    right = (c * this.mNearClippingDistance) / this.mConvergence;

    // Set the Projection Matrix
    this.mProjectionMatrix = m4.orthographic(
      left,
      right,
      bottom,
      top,
      this.mNearClippingDistance,
      this.mFarClippingDistance
    );

    // Displace the world to right
    this.mModelViewMatrix = m4.translation(
      this.mEyeSeparation / 2,
      0.0,
      0.0
    );
  };

  this.ApplyRightFrustum = function() {
    let top, bottom, left, right;
    top = this.mNearClippingDistance * Math.tan(this.mFOV / 2);
    bottom = -top;

    let a = this.mAspectRatio * Math.tan(this.mFOV / 2) * this.mConvergence;
    let b = a - this.mEyeSeparation / 2;
    let c = a + this.mEyeSeparation / 2;

    left = (-c * this.mNearClippingDistance) / this.mConvergence;
    right = (b * this.mNearClippingDistance) / this.mConvergence;

    // Set the Projection Matrix
    this.mProjectionMatrix = m4.orthographic(
      left,
      right,
      bottom,
      top,
      this.mNearClippingDistance,
      this.mFarClippingDistance
    );

    // Displace the world to left
    this.mModelViewMatrix = m4.translation(
      -this.mEyeSeparation / 2,
      0.0,
      0.0
    );
  };
}
let matrix = null;
let vec = {
  alpha: 0,
  beta: 0,
  gamma: 0
}
function requestDeviceOrientation() {
  if (typeof DeviceOrientationEvent !== 'undefined' &&
    typeof DeviceOrientationEvent.requestPermission === 'function') {
    DeviceOrientationEvent.requestPermission()
      .then(response => {
        console.log(response);
        if (response === 'granted') {
          console.log('Permission granted');
          window.addEventListener('deviceorientation', e => {
            matrix = getRotationMatrix(e.alpha, e.beta, e.gamma);
            vec.alpha = e.alpha
            vec.beta = e.beta
            vec.gamma = e.gamma
          }, true);
          animate();

        }
      }).catch((err => {
        console.log('Err', err);
      }));
  } else
    console.log('not iOS');
}

var degtorad = Math.PI / 180; // Degree-to-Radian conversion

function getRotationMatrix(alpha, beta, gamma) {

  var _x = beta ? beta * degtorad : 0; // beta value
  var _y = gamma ? gamma * degtorad : 0; // gamma value
  var _z = alpha ? alpha * degtorad : 0; // alpha value

  var cX = Math.cos(_x);
  var cY = Math.cos(_y);
  var cZ = Math.cos(_z);
  var sX = Math.sin(_x);
  var sY = Math.sin(_y);
  var sZ = Math.sin(_z);

  //
  // ZXY rotation matrix construction.
  //

  var m11 = cZ * cY - sZ * sX * sY;
  var m12 = - cX * sZ;
  var m13 = cY * sZ * sX + cZ * sY;

  var m21 = cY * sZ + cZ * sX * sY;
  var m22 = cZ * cX;
  var m23 = sZ * sY - cZ * cY * sX;

  var m31 = - cX * sY;
  var m32 = sX;
  var m33 = cX * cY;

  return [
    m11, m12, m13, 0,
    m21, m22, m23, 0,
    m31, m32, m33, 0, 0, 0, 0, 1
  ];

};

function vectorFromVec(vector0) {
  // Convert angles to radians
  const alphaRad = (vector0.beta * Math.PI) / 180;
  const betaRad = (vector0.gamma * Math.PI) / 180;
  const gammaRad = (vector0.alpha * Math.PI) / 180;

  // Define the initial vector along the x-axis
  let vector = [0, 0, 1];

  // Rotation around the z-axis (gamma)
  const rotZ = [
    [Math.cos(gammaRad), -Math.sin(gammaRad), 0],
    [Math.sin(gammaRad), Math.cos(gammaRad), 0],
    [0, 0, 1]
  ];
  vector = multiplyMatrixVector(rotZ, vector);

  // Rotation around the y-axis (beta)
  const rotY = [
    [Math.cos(betaRad), 0, Math.sin(betaRad)],
    [0, 1, 0],
    [-Math.sin(betaRad), 0, Math.cos(betaRad)]
  ];
  vector = multiplyMatrixVector(rotY, vector);

  // Rotation around the x-axis (alpha)
  const rotX = [
    [1, 0, 0],
    [0, Math.cos(alphaRad), -Math.sin(alphaRad)],
    [0, Math.sin(alphaRad), Math.cos(alphaRad)]
  ];
  vector = multiplyMatrixVector(rotX, vector);

  return vector;
}

function multiplyMatrixVector(matrix, vector) {
  const result = [];
  for (let i = 0; i < matrix.length; i++) {
    let sum = 0;
    for (let j = 0; j < vector.length; j++) {
      sum += matrix[i][j] * vector[j];
    }
    result.push(sum);
  }
  return result;
}

let context, audio, source, filter, panner;

function playSomeMusic() {
  audio = document.getElementById('someMusic');
  let checkBox = document.getElementById('FLTR')

  audio.addEventListener('play', () => {
    if (!context) {
      context = new AudioContext();
      source = context.createMediaElementSource(audio);
      panner = context.createPanner();
      filter = context.createBiquadFilter();

      source.connect(panner);
      panner.connect(filter);
      filter.connect(context.destination);

      filter.type = 'peaking';
      filter.Q.value = 1;
      filter.frequency.value = 500;
      filter.gain.value = 20;
      context.resume();
    }
  })


  audio.addEventListener('pause', () => {
    console.log('pause');
    context.resume();
  })
  checkBox.addEventListener('change', function() {
    if (checkBox.checked) {
      panner.disconnect();
      panner.connect(filter);
      filter.connect(context.destination);
    } else {
      panner.disconnect();
      panner.connect(context.destination);
    }
  });
}