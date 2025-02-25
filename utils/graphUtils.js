"use strict";
// FMG utils related to graph

// check if new grid graph should be generated or we can use the existing one
function shouldRegenerateGrid(grid, expectedSeed) {
    if (expectedSeed && expectedSeed !== grid.seed) {
      console.log("Seed 不匹配，提前返回 true");
      return true;
    }

    const cellsDesired = +byId("pointsInput").dataset.cells;
    if (cellsDesired !== grid.cellsDesired) {
      console.log("cellsDesired 不匹配，提前返回 true");
      return true;
    }

    const newSpacing = rn(Math.sqrt((graphWidth * graphHeight) / cellsDesired), 2);
    const newCellsX = Math.floor((graphWidth + 0.5 * newSpacing - 1e-10) / newSpacing);
    const newCellsY = Math.floor((graphHeight + 0.5 * newSpacing - 1e-10) / newSpacing);

    console.log("shouldRegenerateGrid: ", grid.spacing !== newSpacing || grid.cellsX !== newCellsX || grid.cellsY !== newCellsY);
    return grid.spacing !== newSpacing || grid.cellsX !== newCellsX || grid.cellsY !== newCellsY;
}



function buildCellNeighbors(cellsV) {
    /**
     * 构建每个 cell 的邻居数组（grid.cells.c）
     * @param {Array} cellsV - 每个 cell 的顶点索引数组，例如 grid.cells.v
     * @returns {Array} neighbors - 邻居数组，每个元素是一个 cell 索引数组
     */
    const vertexToCells = new Map();
    // 1. 遍历每个 cell 的顶点，将顶点映射到 cell 索引
    cellsV.forEach((cellVertices, cellIndex) => {
        cellVertices.forEach(v => {
            if (!vertexToCells.has(v)) {
                vertexToCells.set(v, []);
            }
            vertexToCells.get(v).push(cellIndex);
        });
    });

    // 2. 对于每个 cell，根据其所有顶点找到共享该顶点的其它 cell
    const neighbors = [];
    cellsV.forEach((cellVertices, cellIndex) => {
        const nbSet = new Set();
        cellVertices.forEach(v => {
            vertexToCells.get(v).forEach(otherCellIndex => {
                if (otherCellIndex !== cellIndex) {
                    nbSet.add(otherCellIndex);
                }
            });
        });
        neighbors[cellIndex] = Array.from(nbSet);
    });

    return neighbors;
}


function generateGrid() {
    Math.random = aleaPRNG(seed); // reset PRNG
    const {spacing, cellsDesired, boundary, points, cellsX, cellsY} = placePoints();
    const {cells, vertices} = calculateVoronoi(points, boundary);
    console.log("before buildCellNeighbors, cells", cells, "cells.v", cells.v);
    // cells.c = buildCellNeighbors(cells.v);
    console.log("after buildCellNeighbors, cells", cells, "cells.v:", cells.v);
    return {spacing, cellsDesired, boundary, points, cellsX, cellsY, cells, vertices, seed};
}


function generateHexGrid(width, height, size) {
    const points = [];
    const xSpacing = Math.sqrt(3) * size;  // 水平间距
    const ySpacing = 1.5 * size;           // 垂直间距

    // 计算大致的行数和列数
    const rows = Math.ceil(height / ySpacing);
    const cols = Math.ceil(width / xSpacing);

    for (let row = 0; row < rows; row++) {
        // 对于奇数行添加水平偏移
        const xOffset = (row % 2) * (xSpacing / 2);
        for (let col = 0; col < cols; col++) {
            const x = col * xSpacing + xOffset;
            const y = row * ySpacing;
            // 只收录在画布范围内的点
            if (x < width && y < height) {
                points.push([x, y]);
            }
        }
    }
    return points;
}


// place random points to calculate Voronoi diagram
function placePoints_old() {
    TIME && console.time("placePoints");
    const cellsDesired = +byId("pointsInput").dataset.cells;
    const spacing = rn(Math.sqrt((graphWidth * graphHeight) / cellsDesired), 2); // spacing between points before jirrering

    const boundary = getBoundaryPoints(graphWidth, graphHeight, spacing);
    const points = getJitteredGrid(graphWidth, graphHeight, spacing); // points of jittered square grid


    const cellsX = Math.floor((graphWidth + 0.5 * spacing - 1e-10) / spacing);
    const cellsY = Math.floor((graphHeight + 0.5 * spacing - 1e-10) / spacing);
    TIME && console.timeEnd("placePoints");

    return {spacing, cellsDesired, boundary, points, cellsX, cellsY};
}

function placePoints() {
    TIME && console.time("placePoints");
    const cellsDesired = +byId("pointsInput").dataset.cells;
    // 根据期望 cell 数量和总面积计算一个合适的 size（你可以根据实际需求调整）
    const size = rn(Math.sqrt((graphWidth * graphHeight) / cellsDesired), 2);

    // 生成边界点（如果不需要可以适当调整或保留原有方式）
    const boundary = getBoundaryPoints(graphWidth, graphHeight, size);

    // 用新的六边形网格函数生成点阵
    const points = generateHexGrid(graphWidth, graphHeight, size);

    // 根据六边形网格计算列和行数（这里的计算需要和 generateHexGrid 保持一致）
    const cellsX = Math.ceil(graphWidth / (1.5 * size));
    const cellsY = Math.ceil(graphHeight / (Math.sqrt(3) * size));

    TIME && console.timeEnd("placePoints");
    return {spacing: size, cellsDesired, boundary, points, cellsX, cellsY};
}


// calculate Delaunay and then Voronoi diagram
function calculateVoronoi(points, boundary) {
    TIME && console.time("calculateDelaunay");
    const allPoints = points.concat(boundary);
    const delaunay = Delaunator.from(allPoints);
    TIME && console.timeEnd("calculateDelaunay");

    TIME && console.time("calculateVoronoi");
    const voronoi = new Voronoi(delaunay, allPoints, points.length);

    const cells = voronoi.cells;
    cells.i = createTypedArray({maxValue: points.length, length: points.length}).map((_, i) => i); // array of indexes
    const vertices = voronoi.vertices;
    TIME && console.timeEnd("calculateVoronoi");

    return {cells, vertices};
}

// add points along map edge to pseudo-clip voronoi cells
function getBoundaryPoints(width, height, spacing) {
    const offset = rn(-1 * spacing);
    const bSpacing = spacing * 2;
    const w = width - offset * 2;
    const h = height - offset * 2;
    const numberX = Math.ceil(w / bSpacing) - 1;
    const numberY = Math.ceil(h / bSpacing) - 1;
    const points = [];

    for (let i = 0.5; i < numberX; i++) {
        let x = Math.ceil((w * i) / numberX + offset);
        points.push([x, offset], [x, h + offset]);
    }

    for (let i = 0.5; i < numberY; i++) {
        let y = Math.ceil((h * i) / numberY + offset);
        points.push([offset, y], [w + offset, y]);
    }

    return points;
}

// get points on a regular square grid and jitter them a bit
function getJitteredGrid(width, height, spacing) {
    const radius = spacing / 2; // square radius
    const jittering = radius * 0.9; // max deviation
    const doubleJittering = jittering * 2;
    const jitter = () => Math.random() * doubleJittering - jittering;

    let points = [];
    for (let y = radius; y < height; y += spacing) {
        for (let x = radius; x < width; x += spacing) {
            const xj = Math.min(rn(x + jitter(), 2), width);
            const yj = Math.min(rn(y + jitter(), 2), height);
            points.push([xj, yj]);
        }
    }
    return points;
}


// return cell index on a regular square grid
function findGridCell_old(x, y, grid) {
    return (
        Math.floor(Math.min(y / grid.spacing, grid.cellsY - 1)) * grid.cellsX +
        Math.floor(Math.min(x / grid.spacing, grid.cellsX - 1))
    );
}

function findGridCell_new1(x, y, grid) {
    const col = Math.floor(Math.min(x / grid.spacing, grid.cellsX - 1));
    const row = Math.floor(Math.min(y / grid.spacing, grid.cellsY - 1));
    return row * grid.cellsX + col;
}

// 新的定位函数：返回六边形网格中的 cell 索引
function findGridCell(x, y, grid) {
    const size = grid.spacing;             // 六边形基本尺寸
    const colWidth = 1.5 * size;             // 横向单元格宽度
    const rowHeight = Math.sqrt(3) * size;   // 纵向单元格高度

    // 根据 y 坐标计算行号
    const row = Math.floor(y / rowHeight);
    // 根据行号决定水平偏移（交错排列：奇数行/偶数行不同）
    const xOffset = (row % 2) * (colWidth / 2);
    // 根据 x 坐标（减去偏移量）计算列号
    const col = Math.floor((x - xOffset) / colWidth);

    // 将二维索引转换为一维索引（假设行优先）
    return row * grid.cellsX + col;
}


// return array of cell indexes in radius on a regular square grid
function findGridAll_(x, y, radius) {
    const c = grid.cells.c;
    let r = Math.floor(radius / grid.spacing);
    let found = [findGridCell(x, y, grid)];
    if (!r || radius === 1) return found;
    if (r > 0) found = found.concat(c[found[0]]);
    if (r > 1) {
        let frontier = c[found[0]];
        while (r > 1) {
            let cycle = frontier.slice();
            frontier = [];
            cycle.forEach(function (s) {
                c[s].forEach(function (e) {
                    if (found.indexOf(e) !== -1) return;
                    found.push(e);
                    frontier.push(e);
                });
            });
            r--;
        }
    }

    return found;
}



// 新的 findGridAll：基于 cell 中心点与鼠标位置的距离
function findGridAll(x, y, radius) {
    const selected = [];
    // 遍历 grid.points 数组，假设每个元素是 [centerX, centerY]
    grid.points.forEach((pt, i) => {
        const dx = pt[0] - x;
        const dy = pt[1] - y;
        if (dx * dx + dy * dy <= radius * radius) {
            selected.push(i);
        }
    });
    return selected;
}


// return closest pack points quadtree datum
function find(x, y, radius = Infinity) {
    return pack.cells.q.find(x, y, radius);
}

// return closest cell index
function findCell(x, y, radius = Infinity) {
    if (!pack.cells?.q) return;
    const found = pack.cells.q.find(x, y, radius);
    return found ? found[2] : undefined;
}

// return array of cell indexes in radius
function findAll(x, y, radius) {
    const found = pack.cells.q.findAll(x, y, radius);
    return found.map(r => r[2]);
}

// get polygon points for packed cells knowing cell id
function getPackPolygon(i) {
    return pack.cells.v[i].map(v => pack.vertices.p[v]);
}

// get polygon points for initial cells knowing cell id
function getGridPolygon(i) {
    return grid.cells.v[i].map(v => grid.vertices.p[v]);
}

// mbostock's poissonDiscSampler
function* poissonDiscSampler(x0, y0, x1, y1, r, k = 3) {
    if (!(x1 >= x0) || !(y1 >= y0) || !(r > 0)) throw new Error();

    const width = x1 - x0;
    const height = y1 - y0;
    const r2 = r * r;
    const r2_3 = 3 * r2;
    const cellSize = r * Math.SQRT1_2;
    const gridWidth = Math.ceil(width / cellSize);
    const gridHeight = Math.ceil(height / cellSize);
    const grid = new Array(gridWidth * gridHeight);
    const queue = [];

    function far(x, y) {
        const i = (x / cellSize) | 0;
        const j = (y / cellSize) | 0;
        const i0 = Math.max(i - 2, 0);
        const j0 = Math.max(j - 2, 0);
        const i1 = Math.min(i + 3, gridWidth);
        const j1 = Math.min(j + 3, gridHeight);
        for (let j = j0; j < j1; ++j) {
            const o = j * gridWidth;
            for (let i = i0; i < i1; ++i) {
                const s = grid[o + i];
                if (s) {
                    const dx = s[0] - x;
                    const dy = s[1] - y;
                    if (dx * dx + dy * dy < r2) return false;
                }
            }
        }
        return true;
    }

    function sample(x, y) {
        queue.push((grid[gridWidth * ((y / cellSize) | 0) + ((x / cellSize) | 0)] = [x, y]));
        return [x + x0, y + y0];
    }

    yield sample(width / 2, height / 2);

    pick: while (queue.length) {
        const i = (Math.random() * queue.length) | 0;
        const parent = queue[i];

        for (let j = 0; j < k; ++j) {
            const a = 2 * Math.PI * Math.random();
            const r = Math.sqrt(Math.random() * r2_3 + r2);
            const x = parent[0] + r * Math.cos(a);
            const y = parent[1] + r * Math.sin(a);
            if (0 <= x && x < width && 0 <= y && y < height && far(x, y)) {
                yield sample(x, y);
                continue pick;
            }
        }

        const r = queue.pop();
        if (i < queue.length) queue[i] = r;
    }
}

// filter land cells
function isLand(i) {
    return pack.cells.h[i] >= 20;
}

// filter water cells
function isWater(i) {
    return pack.cells.h[i] < 20;
}

// findAll d3.quandtree search from https://bl.ocks.org/lwthatcher/b41479725e0ff2277c7ac90df2de2b5e
void (function addFindAll() {
    const Quad = function (node, x0, y0, x1, y1) {
        this.node = node;
        this.x0 = x0;
        this.y0 = y0;
        this.x1 = x1;
        this.y1 = y1;
    };

    const tree_filter = function (x, y, radius) {
        const t = {x, y, x0: this._x0, y0: this._y0, x3: this._x1, y3: this._y1, quads: [], node: this._root};
        if (t.node) t.quads.push(new Quad(t.node, t.x0, t.y0, t.x3, t.y3));
        radiusSearchInit(t, radius);

        var i = 0;
        while ((t.q = t.quads.pop())) {
            i++;

            // Stop searching if this quadrant can’t contain a closer node.
            if (
                !(t.node = t.q.node) ||
                (t.x1 = t.q.x0) > t.x3 ||
                (t.y1 = t.q.y0) > t.y3 ||
                (t.x2 = t.q.x1) < t.x0 ||
                (t.y2 = t.q.y1) < t.y0
            )
                continue;

            // Bisect the current quadrant.
            if (t.node.length) {
                t.node.explored = true;
                var xm = (t.x1 + t.x2) / 2,
                    ym = (t.y1 + t.y2) / 2;

                t.quads.push(
                    new Quad(t.node[3], xm, ym, t.x2, t.y2),
                    new Quad(t.node[2], t.x1, ym, xm, t.y2),
                    new Quad(t.node[1], xm, t.y1, t.x2, ym),
                    new Quad(t.node[0], t.x1, t.y1, xm, ym)
                );

                // Visit the closest quadrant first.
                if ((t.i = ((y >= ym) << 1) | (x >= xm))) {
                    t.q = t.quads[t.quads.length - 1];
                    t.quads[t.quads.length - 1] = t.quads[t.quads.length - 1 - t.i];
                    t.quads[t.quads.length - 1 - t.i] = t.q;
                }
            }

            // Visit this point. (Visiting coincident points isn’t necessary!)
            else {
                var dx = x - +this._x.call(null, t.node.data),
                    dy = y - +this._y.call(null, t.node.data),
                    d2 = dx * dx + dy * dy;
                radiusSearchVisit(t, d2);
            }
        }
        return t.result;
    };
    d3.quadtree.prototype.findAll = tree_filter;

    var radiusSearchInit = function (t, radius) {
        t.result = [];
        (t.x0 = t.x - radius), (t.y0 = t.y - radius);
        (t.x3 = t.x + radius), (t.y3 = t.y + radius);
        t.radius = radius * radius;
    };

    var radiusSearchVisit = function (t, d2) {
        t.node.data.scanned = true;
        if (d2 < t.radius) {
            do {
                t.result.push(t.node.data);
                t.node.data.selected = true;
            } while ((t.node = t.node.next));
        }
    };
})();

// draw raster heightmap preview (not used in main generation)
function drawHeights({heights, width, height, scheme, renderOcean}) {
    const canvas = document.createElement("canvas");
    canvas.width = width;
    canvas.height = height;
    const ctx = canvas.getContext("2d");
    const imageData = ctx.createImageData(width, height);

    const getHeight = height => (height < 20 ? (renderOcean ? height : 0) : height);

    for (let i = 0; i < heights.length; i++) {
        const color = scheme(1 - getHeight(heights[i]) / 100);
        const {r, g, b} = d3.color(color);

        const n = i * 4;
        imageData.data[n] = r;
        imageData.data[n + 1] = g;
        imageData.data[n + 2] = b;
        imageData.data[n + 3] = 255;
    }

    ctx.putImageData(imageData, 0, 0);
    return canvas.toDataURL("image/png");
}
