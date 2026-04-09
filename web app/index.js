
let zoom = 1.00;
let zMin = 0.05;
let zMax = 10.00;
let sensitivity = 0.001;

let offsetX = 0;
let offsetY = 0;

let particles = [];

class particle {
    constructor(mass, radius, x, y, Vx, Vy) {
        this.mass = mass;
        this.radius = radius;
        this.x = x;
        this.y = y;
        // cpp constuctor here
        particles.push(this);
    }

    update() {
        // cpp retreiver here
        this.x = 0;
        this.y = 0;
    }

    draw() {
        circle(this.x, this.y, this.radius);
    }
}

document.addEventListener('contextmenu', function (event) {
    event.preventDefault();
});

function setup() {
    createCanvas(windowWidth, windowHeight);
    offsetX = width / 2;
    offsetY = height / 2;
}


function draw() {
    background(220);

    push();
    translate(offsetX, offsetY);
    scale(zoom);

    drawGrid();
    fill(255, 100, 100);
    rectMode(CENTER);

    particles.forEach(part => {
        part.draw();
    });

    pop();
}

function mouseWheel(event) {
    let s = 1 - event.delta * sensitivity;
    let nextZoom = constrain(zoom * s, zMin, zMax);

    offsetX -= (mouseX - offsetX) * (nextZoom / zoom - 1);
    offsetY -= (mouseY - offsetY) * (nextZoom / zoom - 1);

    zoom = nextZoom;

    return false;
}

function mouseDragged() {
    if (mouseButton === RIGHT) {
        offsetX += mouseX - pmouseX;
        offsetY += mouseY - pmouseY;
    }
}
let mx = 0;
let my = 0;

function mousePressed() {
    if (mouseButton === LEFT) {
        mx = mouseX;
        my = mouseY;
    }
}

function mouseReleased() {
    if (mouseButton === LEFT) {
        let Vx = mouseX - mx;
        let Vy = mouseY - my;
        new particle(100, 100, (mx - offsetX) / zoom, (my - offsetY) / zoom, Vx, Vy);
    }
}

function drawGrid() {
    stroke(200);
    strokeWeight(1 / zoom);

    let step = 50;

    let nw = -offsetX / zoom;
    let se = (width - offsetX) / zoom;
    let ne = -offsetY / zoom;
    let sw = (height - offsetY) / zoom;

    for (let x = Math.floor(nw / step) * step; x < se; x += step) {
        line(x, ne, x, sw);
    }

    for (let y = Math.floor(ne / step) * step; y < sw; y += step) {
        line(nw, y, se, y);
    }
}
function isVisible(x, y, w, h) {
    let left = -offsetX / zoom;
    let right = (width - offsetX) / zoom;
    let top = -offsetY / zoom;
    let bottom = (height - offsetY) / zoom;

    return (x + w > left && x < right && y + h > top && y < bottom);
}

function windowResized() {
    resizeCanvas(windowWidth, windowHeight);
}