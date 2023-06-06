'use strict';

let gl;                         // The webgl context.
let surface;                    // A surface model
let shProgram;                  // A shader program
let spaceball;                  // A SimpleRotator object that lets the user rotate the view by mouse.
let sphere;
let userPointCoord;
let userScaleFactor;
let video;
let texture;
let texturevid;
let track;
let webCam;

function deg2rad(angle) {
    return angle * Math.PI / 180;
}


// Constructor
function Model(name) {
    this.name = name;
    this.iVertexBuffer = gl.createBuffer();
    this.iTextureBuffer = gl.createBuffer();
    this.count = 0;
    this.countT = 0;

    this.BufferData = function (vertices) {

        gl.bindBuffer(gl.ARRAY_BUFFER, this.iVertexBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STREAM_DRAW);

        this.count = vertices.length / 3;
    }

    this.TextureBufferData = function (points) {

        gl.bindBuffer(gl.ARRAY_BUFFER, this.iTextureBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(points), gl.STREAM_DRAW);

        this.countT = points.length / 2;
    }

    this.Draw = function () {

        gl.bindBuffer(gl.ARRAY_BUFFER, this.iVertexBuffer);
        gl.vertexAttribPointer(shProgram.iAttribVertex, 3, gl.FLOAT, false, 0, 0);
        gl.enableVertexAttribArray(shProgram.iAttribVertex);

        gl.bindBuffer(gl.ARRAY_BUFFER, this.iTextureBuffer);
        gl.vertexAttribPointer(shProgram.iAttribTexture, 2, gl.FLOAT, false, 0, 0);
        gl.enableVertexAttribArray(shProgram.iAttribTexture);

        gl.drawArrays(gl.TRIANGLE_STRIP, 0, this.count);
    }

    this.DrawSphere = function () {
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

    this.iUserPoint = -1;
    this.iScale = 1.0;
    this.iUP = -1;

    this.Use = function () {
        gl.useProgram(this.prog);
    }
}

function getWebcam() {
    navigator.getUserMedia({ video: true, audio: false }, function (stream) {
        video.srcObject = stream;
        track = stream.getTracks()[0];
    }, function (e) {
        console.error('Rejected!', e);
    });
}

function CreateWebCamTexture() {
    texturevid = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texturevid);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
}


/* Draws a colored cube, along with a set of coordinate axes.
 * (Note that the use of the above drawPrimitive function is not an efficient
 * way to draw with WebGL.  Here, the geometry is so simple that it doesn't matter.)
 */
function draw() {
    let D = document;
    let spans = D.getElementsByClassName("slider-value");

    gl.clearColor(0, 0, 0, 1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    /* Set the values of the projection transformation */
    let projection = m4.perspective(Math.PI / 8, 1, 8, 12);
    let conv, // convergence
        eyes, // eye separation
        ratio, // aspect ratio
        fov; // field of view
    conv = 1500.0;
    conv = D.getElementById("conv").value;
    spans[3].innerHTML=conv;
    eyes = 12.0;
    eyes = D.getElementById("eyes").value;
    spans[0].innerHTML=eyes;
    ratio = 1.0;
    fov = 0.45;
    fov = D.getElementById("fov").value;
    spans[1].innerHTML=fov;
    let top, bottom, left, right, near, far;
    near = 30.0;
    near = D.getElementById("near").value-0.0;
    spans[2].innerHTML=near;
    far = 2000.0;

    top = near * Math.tan(fov / 2.0);
    bottom = -top;

    let a = ratio * Math.tan(fov / 2.0) * conv;

    let b = a - eyes / 2;
    let c = a + eyes / 2;

    left = -b * near / conv;
    right = c * near / conv;

    // console.log(left, right, bottom, top, near, far);

    let projectionLeft = m4.orthographic(left, right, bottom, top, near, far);

    left = -c * near / conv;
    right = b * near / conv;

    let projectionRight = m4.orthographic(left, right, bottom, top, near, far);

    /* Get the view matrix from the SimpleRotator object.*/
    let modelView = spaceball.getViewMatrix();
    let initialMat = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]

    let rotateToPointZero = m4.axisRotation([0.707, 0.707, 0], 0.);
    let translateToPointZero = m4.translation(0.0, 0, -20);
    let translateToLeft = m4.translation(-0.03, 0, -25);
    let translateToRight = m4.translation(0.03, 0, -25);
    let translateToBack = m4.translation(-5.0, -5, -25);

    let matAccum = m4.multiply(rotateToPointZero, modelView);
    let matAccum0 = m4.multiply(rotateToPointZero, initialMat);
	if (matrix != null) {
		matAccum = m4.multiply(rotateToPointZero, matrix);
	}
    let matAccum1 = m4.multiply(translateToPointZero, matAccum);
    let matAccum01 = m4.multiply(translateToPointZero, matAccum0);
    let matAccumLeft = m4.multiply(translateToLeft, matAccum1);
    let matAccumRight = m4.multiply(translateToRight, matAccum1);
    let matAccumBack = m4.multiply(translateToBack, m4.multiply(matAccum01,m4.scaling(10,10,10)));

    /* Multiply the projection matrix times the modelview matrix to give the
       combined transformation matrix, and send that to the shader program. */
    let modelViewProjection = m4.multiply(projection, matAccum1);
    let modelViewProjection0 = m4.multiply(projection, matAccum01);
    
    
    gl.uniformMatrix4fv(shProgram.iModelViewMatrix, false, matAccumBack);
    gl.uniformMatrix4fv(shProgram.iProjectionMatrix, false, [0.14236292243003845, 0, 0, 0, 0, 0.14236292243003845, 0, 0, 0, 0, -0.0010152284521609545, 0, -0.01708354987204075, -0, -1.0304569005966187, 1]);
    console.log(projectionLeft)


    gl.uniform1f(shProgram.iScale, 1)
	gl.enable(gl.TEXTURE_2D);
    /*gl.bindTexture(gl.TEXTURE_2D, texturevid);
    gl.texImage2D(
        gl.TEXTURE_2D,
        0,
        gl.RGBA,
        gl.RGBA,
        gl.UNSIGNED_BYTE,
        video
    );
    webCam.Draw();*/

    gl.clear(gl.DEPTH_BUFFER_BIT);

    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.uniformMatrix4fv(shProgram.iModelViewProjectionMatrix, false, modelViewProjection);
    gl.uniform1i(shProgram.iTMU, 0);
    gl.enable(gl.TEXTURE_2D);
    gl.uniform1f(shProgram.iScale, userScaleFactor)
    gl.uniformMatrix4fv(shProgram.iModelViewMatrix, false, matAccumLeft);
    gl.uniformMatrix4fv(shProgram.iProjectionMatrix, false, projectionLeft);
    gl.colorMask(true, false, false, false);
    surface.Draw();

    gl.clear(gl.DEPTH_BUFFER_BIT);

    gl.uniformMatrix4fv(shProgram.iModelViewMatrix, false, matAccumRight);
    gl.uniformMatrix4fv(shProgram.iProjectionMatrix, false, projectionRight);
    gl.colorMask(false, true, true, false);
    surface.Draw();

    gl.clear(gl.DEPTH_BUFFER_BIT);

    gl.colorMask(true, true, true, true);


    
}

function loop(){
    draw()
    window.requestAnimationFrame(loop)

}

function conjugation(z, b, a, c) {
    let r = a * (1 - Math.cos(Math.PI * 2 * z / c)) + 0.8 * 1.5;
    let x = 0.75 * r * Math.cos(b);
    let y = 0.75 * r * Math.sin(b)
    let z1 = 0.75 * z
    return { x: x, y: y, z: z1 }
}

function CreateSurfaceData() {
    let vertexList = [];
    let i = 0;
    let j = 0;
    let a = -0.2633257397764612;
    let c = 3.2;
    let b = 1.6099263856487789;
    let step = 0.05;

    while (i < b + 0.5) {
        while (j < Math.PI * 2) {
            let v1 = conjugation(i, j, a, c)
            let v2 = conjugation(i + step, j, a, c)
            let v3 = conjugation(i, j + step, a, c)
            let v4 = conjugation(i + step, j + step, a, c)
            vertexList.push(v1.x, v1.y, v1.z);
            vertexList.push(v2.x, v2.y, v2.z);
            vertexList.push(v3.x, v3.y, v3.z);
            vertexList.push(v2.x, v2.y, v2.z);
            vertexList.push(v4.x, v4.y, v4.z);
            vertexList.push(v3.x, v3.y, v3.z);
            j += step
        }
        j = 0;
        i += step
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
    let step = 0.05;

    while (i < b + 0.5) {
        while (j < Math.PI * 2) {
            let u = map(i, 0, b, 0, 1);
            let v = map(j, 0, Math.PI * 2, 0, 1);
            texCoordList.push(u, v);
            u = map(i + step, 0, b, 0, 1);
            texCoordList.push(u, v);
            u = map(i, 0, b, 0, 1);
            v = map(j + step, 0, Math.PI * 2, 0, 1);
            texCoordList.push(u, v);
            u = map(i + step, 0, b, 0, 1);
            v = map(j, 0, Math.PI * 2, 0, 1);
            texCoordList.push(u, v);
            u = map(i + step, 0, Math.PI, 0, 1);
            v = map(j + step, 0, Math.PI * 2, 0, 1);
            texCoordList.push(u, v);
            u = map(i, 0, b, 0, 1);
            v = map(j + step, 0, Math.PI * 2, 0, 1);
            texCoordList.push(u, v);
            j += step;
        }
        j = 0
        i += step;
    }
    return texCoordList;
}


function cassini(u, z1) {
    let x = 0.6 * Math.sqrt(r(u, z1)) * Math.cos(u)
    let y = 0.6 * Math.sqrt(r(u, z1)) * Math.sin(u)
    let z = 0.6 * z1;
    return { x: x, y: y, z: z }
}

function r(u, z) {
    let a = 1.0;
    let b = a + z;
    return (a ** 2 * Math.cos(2 * u) + Math.sqrt((b ** 4 - a ** 4) + a ** 4 * Math.cos(2 * u) ** 2));
}

function CreateSphereSurface(r = 0.05) {
    let vertexList = [];
    let lon = -Math.PI;
    let lat = -Math.PI * 0.5;
    while (lon < Math.PI) {
        while (lat < Math.PI * 0.5) {
            let v1 = sphereSurfaceData(r, lon, lat);
            vertexList.push(v1.x, v1.y, v1.z);
            lat += 0.05;
        }
        lat = -Math.PI * 0.5
        lon += 0.05;
    }
    return vertexList;
}

function sphereSurfaceData(r, u, v) {
    let x = r * Math.sin(u) * Math.cos(v);
    let y = r * Math.sin(u) * Math.sin(v);
    let z = r * Math.cos(u);
    return { x: x, y: y, z: z };
}

/* Initialize the WebGL context. Called from init() */
function initGL() {
    let prog = createProgram(gl, vertexShaderSource, fragmentShaderSource);

    shProgram = new ShaderProgram('Basic', prog);
    shProgram.Use();

    shProgram.iAttribVertex = gl.getAttribLocation(prog, "vertex");
    shProgram.iAttribTexture = gl.getAttribLocation(prog, "texture");
    shProgram.iModelViewProjectionMatrix = gl.getUniformLocation(prog, "ModelViewProjectionMatrix");
    shProgram.iModelViewMatrix = gl.getUniformLocation(prog, "ModelViewMatrix");
    shProgram.iProjectionMatrix = gl.getUniformLocation(prog, "ProjectionMatrix");
    shProgram.iTMU = gl.getUniformLocation(prog, 'tmu');
    shProgram.iUserPoint = gl.getUniformLocation(prog, 'userPoint');
    shProgram.iScale = gl.getUniformLocation(prog, 'scl');
    shProgram.iUP = gl.getUniformLocation(prog, 'translateUP');

    surface = new Model('Surface');
    surface.BufferData(CreateSurfaceData());
    LoadTexture();
    surface.TextureBufferData(CreateTextureData());
    sphere = new Model('Sphere');
    sphere.BufferData(CreateSphereSurface())
    webCam = new Model('WC');
    let s =1;
    webCam.BufferData([0.0, 0.0, 0.0, s, 0.0, 0.0, s, s, 0.0, s, s, 0.0, 0.0, s, 0.0, 0.0, 0.0, 0.0]);
    webCam.TextureBufferData([ 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]);

    gl.enable(gl.DEPTH_TEST);
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
    userPointCoord = { x: 0.5, y: 0.5 }
    userScaleFactor = 1.0;
    let canvas;
    try {
        let resolution = Math.min(window.innerHeight, window.innerWidth);
        canvas = document.querySelector('canvas');
        gl = canvas.getContext("webgl");
        canvas.width = resolution;
        canvas.height = resolution;
        gl.viewport(0, 0, resolution, resolution);
        //video = document.createElement('video');
        //video.setAttribute('autoplay', true);
        //window.vid = video;
        //getWebcam();
        //CreateWebCamTexture();
        if (!gl) {
            throw "Browser does not support WebGL";
        }
    }
    catch (e) {
        document.querySelector('"canvas-holder"').innerHTML =
            "<p>Sorry, could not get a WebGL graphics context.</p>";
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

    loop();
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
    texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);

    // gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 512, 512, 0, );

    const image = new Image();
    image.crossOrigin = 'anonymus';
    image.src =  "https://raw.githubusercontent.com/anton-taranec/GW/main/texture.jpg";
    image.onload = () => {
        gl.bindTexture(gl.TEXTURE_2D, texture);
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

let matrix = null;
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
          }, true);
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
