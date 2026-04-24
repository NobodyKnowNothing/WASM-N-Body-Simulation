

function initGraph1() {
    let trace1 = {
        x: [],
        y: [],
        mode: 'lines',
        line: { color: 'red', width: 2 },
        name: 'Mean',
        yaxis: 'y'
    };

    let trace2 = {
        x: [],
        y: [],
        mode: 'lines',
        line: { color: 'blue', width: 2 },
        name: 'Variance',
        yaxis: 'y2'
    };

    let trace3 = {
        x: [],
        y: [],
        mode: 'lines',
        line: { color: 'green', width: 2 },
        name: 'Std. Dev.',
        yaxis: 'y3'
    };

    let layout = {
        legend: {
            x: 1.1
        },
        title: { text: 'Velocity Statistics (m/s)' },
        xaxis: { title: { text: 'Frames rendered' } },
        yaxis: {
            title: { text: 'Mean' },
            autorange: true // Let it scale automatically at first
        },
        yaxis2: {
            title: { text: 'Variance' },
            overlaying: 'y',
            side: 'top',
            anchor: 'free',
            position: 0.5,
        },
        yaxis3: {
            title: { text: 'Standard Deviation' },
            overlaying: 'y',
            side: 'right'
        },
        hovermode: false,
        margin: { t: 40, b: 60, l: 60, r: 30 }
    };

    Plotly.newPlot('chart-container1', [trace1, trace2, trace3], layout, { responsive: true, displayModeBar: false });
}

const sketch = (p) => {
    let zoom = 1.00;
    let zMin = 0.05;
    let zMax = 10.00;
    let sensitivity = 0.001;



    let offsetX = 0;
    let offsetY = 0;

    let dt = 0.1

    let particles = [];

    let statsDiv = null;
    let statsButton = null;
    let statsPlotBool = false;

    let canvas = null;
    class particle {
        constructor(mass, radius, x, y, Vx, Vy) {
            this.mass = mass;
            this.radius = radius;
            this.x = x;
            this.y = y;
            Module._add_particle_(x, y, Vx, Vy, mass, radius);
            particles.push(this);
            this.index = (particles.length - 1);
            Module._setup_verlet_(this.index);
        }

        update() {
            this.x = Module._particle_get_x_(this.index);
            this.y = Module._particle_get_y_(this.index);
        }

        draw() {
            p.fill(255, 100, 100);
            p.circle(this.x, this.y, this.radius * 2);
        }
    }

    document.addEventListener('contextmenu', function (event) {
        event.preventDefault();
    });

    p.statsPlot = () => {
        if (!statsPlotBool) {
            statsDiv = p.createDiv('');
            statsDiv.position(p.windowWidth * 0.03, p.windowHeight * 0.06);
            statsDiv.id('chart-container1')
            statsDiv.size(p.windowWidth * 0.4, p.windowHeight * 0.3);
            initGraph1();
            statsPlotBool = true;
            statsButton.style('border-style', 'inset');
            statsButton.style('background-color', '#7a7a7aff');
        } else {
            statsDiv.remove();
            statsDiv = null;
            statsPlotBool = false;
            statsButton.style('border-style', 'outset');
            statsButton.style('background-color', '#ffffffbb');

        }
    }


    p.setup = () => {
        p.createCanvas(p.windowWidth, p.windowHeight);

        statsButton = p.createButton('Stats');
        statsButton.position(p.windowWidth * 0.03, p.windowHeight * 0.03);
        statsButton.mousePressed(p.statsPlot);
        statsButton.style('appearance', 'none');
        statsButton.style('-webkit-appearance', 'none'); // Needed for Safari/Chrome

        offsetX = p.width / 2;
        offsetY = p.height / 2;

    }

    p.draw = () => {
        if (!Module || !Module._add_particle_) return;
        p.background(220);

        p.push();
        p.translate(offsetX, offsetY);
        p.scale(zoom);

        p.drawGrid();
        p.fill(255, 100, 100);
        p.rectMode(p.CENTER);

        Module._verlet_(dt);
        particles.forEach(part => {
            part.draw();
            part.update();
        });
        if (document.getElementById('chart-container1')) {
            console.log(Module._mean_vel_(), Module._variance_vel_(), Module._std_dev_vel_());
            Plotly.extendTraces('chart-container1', {
                y: [[Module._mean_vel_()], [Module._variance_vel_()], [Module._std_dev_vel_()]],
                x: [[p.frameCount], [p.frameCount], [p.frameCount]]
            }, [0, 1, 2], 300);
        }

        p.pop();
    }

    p.mouseWheel = (event) => {
        let s = 1 - event.delta * sensitivity;
        let nextZoom = p.constrain(zoom * s, zMin, zMax);

        offsetX -= (p.mouseX - offsetX) * (nextZoom / zoom - 1);
        offsetY -= (p.mouseY - offsetY) * (nextZoom / zoom - 1);

        zoom = nextZoom;

        return false;
    }

    p.mouseDragged = () => {
        if (p.mouseButton === p.RIGHT) {
            offsetX += p.mouseX - p.pmouseX;
            offsetY += p.mouseY - p.pmouseY;
        }
    }
    let mx = 0;
    let my = 0;

    p.mousePressed = (event) => {
        if (event.target.tagName === 'CANVAS' && p.mouseButton === p.LEFT) {
            mx = p.mouseX;
            my = p.mouseY;
        }
    }

    p.mouseReleased = (event) => {
        if (event.target.tagName === 'CANVAS' && p.mouseButton === p.LEFT) {
            let Vx = (mx - p.mouseX) / zoom;
            let Vy = (my - p.mouseY) / zoom;
            new particle(100000000000000, 10, (p.mouseX - offsetX) / zoom, (p.mouseY - offsetY) / zoom, Vx, Vy);
        }
    }

    p.drawGrid = () => {
        p.stroke(200);
        p.strokeWeight(1 / zoom);

        let step = 50;

        let nw = -offsetX / zoom;
        let se = (p.width - offsetX) / zoom;
        let ne = -offsetY / zoom;
        let sw = (p.height - offsetY) / zoom;

        for (let x = Math.floor(nw / step) * step; x < se; x += step) {
            p.line(x, ne, x, sw);
        }

        for (let y = Math.floor(ne / step) * step; y < sw; y += step) {
            p.line(nw, y, se, y);
        }
    }
    p.isVisible = (x, y, w, h) => {
        let left = -offsetX / zoom;
        let right = (p.width - offsetX) / zoom;
        let top = -offsetY / zoom;
        let bottom = (p.height - offsetY) / zoom;

        return (x + w > left && x < right && y + h > top && y < bottom);
    }

    p.windowResized = () => {
        p.resizeCanvas(p.windowWidth, p.windowHeight);
    }
}
Module.onRuntimeInitialized = () => {
    new p5(sketch);
}