const Model = (() => {
    const PRIMARIES = {
        R: { x: 0.64, y: 0.33 },
        G: { x: 0.30, y: 0.60 },
        B: { x: 0.15, y: 0.06 }
    };

    function matMul3x3(A, v) {
        return [
            A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2],
            A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2],
            A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2]
        ];
    }

    function matInv3x3(A) {
        const [a,b,c,d,e,f,g,h,i] = [A[0][0],A[0][1],A[0][2], A[1][0],A[1][1],A[1][2], A[2][0],A[2][1],A[2][2]];
        const det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
        if (Math.abs(det) < 1e-12) throw new Error('det zero');
        return [
            [(e*i - f*h)/det, (c*h - b*i)/det, (b*f - c*e)/det],
            [(f*g - d*i)/det, (a*i - c*g)/det, (c*d - a*f)/det],
            [(d*h - e*g)/det, (b*g - a*h)/det, (a*e - b*d)/det]
        ];
    }

    function computeRgbToXyzMatrix(Xw, Yw, Zw) {
        const { R, G, B } = PRIMARIES;
        const M = [
            [ R.x/R.y, G.x/G.y, B.x/B.y ],
            [ 1,       1,       1       ],
            [ (1-R.x-R.y)/R.y, (1-G.x-G.y)/G.y, (1-B.x-B.y)/B.y ]
        ];
        const invM = matInv3x3(M);
        const [Sr, Sg, Sb] = matMul3x3(invM, [Xw, Yw, Zw]);
        return [
            [ M[0][0]*Sr, M[0][1]*Sg, M[0][2]*Sb ],
            [ M[1][0]*Sr, M[1][1]*Sg, M[1][2]*Sb ],
            [ M[2][0]*Sr, M[2][1]*Sg, M[2][2]*Sb ]
        ];
    }

    function getIlluminantData(illuminant) {
        const whitePoints = {
            D65: { X: 95.047, Y: 100, Z: 108.883 },
            D50: { X: 96.421, Y: 100, Z: 82.518 },
            E:   { X: 100,    Y: 100, Z: 100    }
        };
        const w = whitePoints[illuminant];
        if (!w) throw new Error('Unknown illuminant');
        const rgb2xyz = computeRgbToXyzMatrix(w.X, w.Y, w.Z);
        const xyz2rgb = matInv3x3(rgb2xyz);
        return { rgb2xyz, xyz2rgb, white: w };
    }

    function rgbToXyz(r, g, b, rgb2xyz) {
        const lin = (c) => {
            c = c / 255;
            return (c <= 0.04045) ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
        };
        const [R, G, B] = [lin(r), lin(g), lin(b)];
        const [X, Y, Z] = matMul3x3(rgb2xyz, [R, G, B]);
        return { X, Y, Z };
    }

    function xyzToRgb(X, Y, Z, xyz2rgb, strategy) {
        const [Rn, Gn, Bn] = matMul3x3(xyz2rgb, [X, Y, Z]);
        const invLin = (c) => {
            if (c <= 0.0031308) return 12.92 * c;
            return 1.055 * Math.pow(c, 1/2.4) - 0.055;
        };
        let r = invLin(Rn) * 255;
        let g = invLin(Gn) * 255;
        let b = invLin(Bn) * 255;

        let clipped = false;
        if (strategy === 'clip') {
            if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255) clipped = true;
            r = Math.min(255, Math.max(0, r));
            g = Math.min(255, Math.max(0, g));
            b = Math.min(255, Math.max(0, b));
        } else if (strategy === 'scale') {
            const vals = [r, g, b];
            const minV = Math.min(...vals);
            const maxV = Math.max(...vals);
            if (minV < 0 || maxV > 255) {
                clipped = true;
                const range = maxV - minV;
                if (range < 1e-12) {
                    r = g = b = 127.5;
                } else {
                    r = ((r - minV) / range) * 255;
                    g = ((g - minV) / range) * 255;
                    b = ((b - minV) / range) * 255;
                }
            }
        }
        return { r, g, b, clipped };
    }

    function xyzToLab(X, Y, Z, white) {
        const delta = 6 / 29;
        const f = (t) => (t > delta ** 3)
            ? Math.cbrt(t)
            : t / (3 * delta ** 2) + 4 / 29;
        const x = X / white.X;
        const y = Y / white.Y;
        const z = Z / white.Z;
        const fx = f(x), fy = f(y), fz = f(z);
        const L = 116 * fy - 16;
        const a = 500 * (fx - fy);
        const b = 200 * (fy - fz);
        return { L, a, b };
    }

    function labToXyz(L, a, b, white) {
        const fy = (L + 16) / 116;
        const fx = fy + a / 500;
        const fz = fy - b / 200;
        const delta = 6 / 29;
        const fInv = (t) => (t > delta)
            ? t ** 3
            : 3 * delta ** 2 * (t - 4 / 29);
        return {
            X: fInv(fx) * white.X,
            Y: fInv(fy) * white.Y,
            Z: fInv(fz) * white.Z
        };
    }

    function rgbToHsl(r, g, b) {
        const R = r/255, G = g/255, B = b/255;
        const max = Math.max(R,G,B), min = Math.min(R,G,B);
        let h, s, l = (max + min) / 2;
        if (max === min) {
            h = s = 0;
        } else {
            const d = max - min;
            s = (l > 0.5) ? d / (2 - max - min) : d / (max + min);
            switch (max) {
                case R: h = ((G - B) / d + (G < B ? 6 : 0)) / 6; break;
                case G: h = ((B - R) / d + 2) / 6; break;
                case B: h = ((R - G) / d + 4) / 6; break;
            }
        }
        return { h: h * 360, s, l };
    }

    function hslToRgb(h, s, l) {
        const hue = ((h % 360) + 360) % 360 / 360;
        const c = (1 - Math.abs(2*l - 1)) * s;
        const x = c * (1 - Math.abs((hue * 6) % 2 - 1));
        const m = l - c/2;
        let r, g, b;
        if (hue < 1/6) { r=c; g=x; b=0; }
        else if (hue < 2/6) { r=x; g=c; b=0; }
        else if (hue < 3/6) { r=0; g=c; b=x; }
        else if (hue < 4/6) { r=0; g=x; b=c; }
        else if (hue < 5/6) { r=x; g=0; b=c; }
        else { r=c; g=0; b=x; }
        return { r: (r+m)*255, g: (g+m)*255, b: (b+m)*255 };
    }

    function hexToRgb(hex) {
        const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
        return result ? {
            r: parseInt(result[1], 16),
            g: parseInt(result[2], 16),
            b: parseInt(result[3], 16)
        } : null;
    }

    function rgbToHex(r, g, b) {
        const clamp = (v) => Math.round(Math.min(255, Math.max(0, v)));
        return '#' + [clamp(r), clamp(g), clamp(b)].map(c => c.toString(16).padStart(2, '0')).join('');
    }

    return {
        getIlluminantData,
        rgbToXyz,
        xyzToRgb,
        xyzToLab,
        labToXyz,
        rgbToHsl,
        hslToRgb,
        hexToRgb,
        rgbToHex
    };
})();
function runTests() {
    const results = [];
    const assert = (cond, msg) => {
        results.push({ ok: cond, msg });
    };

    const dataD65 = Model.getIlluminantData('D65');
    const xyz1 = Model.rgbToXyz(255,0,0, dataD65.rgb2xyz);
    assert(Math.abs(xyz1.X - 41.2456) < 0.01 && Math.abs(xyz1.Y - 21.2673) < 0.01 && Math.abs(xyz1.Z - 1.9334) < 0.01,
        'RGB(255,0,0) -> XYZ (D65) совпадает с эталоном');
    const lab1 = Model.xyzToLab(xyz1.X, xyz1.Y, xyz1.Z, dataD65.white);
    assert(Math.abs(lab1.L - 53.2408) < 0.02 && Math.abs(lab1.a - 80.0925) < 0.02 && Math.abs(lab1.b - 67.2032) < 0.02,
        'XYZ -> LAB (D65) совпадает с эталоном для красного');
    const rgbBack = Model.xyzToRgb(xyz1.X, xyz1.Y, xyz1.Z, dataD65.xyz2rgb, 'clip');
    assert(Math.abs(rgbBack.r - 255) < 0.5 && Math.abs(rgbBack.g) < 0.5 && Math.abs(rgbBack.b) < 0.5,
        'RGB(255,0,0) -> XYZ -> RGB (D65, clip) даёт исходный');
    const xyz2 = Model.labToXyz(lab1.L, lab1.a, lab1.b, dataD65.white);
    assert(Math.abs(xyz2.X - xyz1.X) < 0.1 && Math.abs(xyz2.Y - xyz1.Y) < 0.1 && Math.abs(xyz2.Z - xyz1.Z) < 0.1,
        'XYZ -> LAB -> XYZ (D65) даёт исходный (погрешность <0.1)');
    const hsl1 = Model.rgbToHsl(255,0,0);
    const rgb3 = Model.hslToRgb(hsl1.h, hsl1.s, hsl1.l);
    assert(Math.abs(rgb3.r - 255) < 0.5 && Math.abs(rgb3.g) < 0.5 && Math.abs(rgb3.b) < 0.5,
        'RGB(255,0,0) -> HSL -> RGB даёт исходный');
 const dataD50 = Model.getIlluminantData('D50');
    const xyzWhite = Model.rgbToXyz(255,255,255, dataD50.rgb2xyz);
    const rgbBackWhite = Model.xyzToRgb(xyzWhite.X, xyzWhite.Y, xyzWhite.Z, dataD50.xyz2rgb, 'clip');
    assert(Math.abs(rgbBackWhite.r - 255) < 0.5 && Math.abs(rgbBackWhite.g - 255) < 0.5 && Math.abs(rgbBackWhite.b - 255) < 0.5,
        'RGB(255,255,255) -> XYZ (D50) -> RGB (D50, clip) даёт исходный');
    const rgbScale = Model.xyzToRgb(100,100,100, dataD65.xyz2rgb, 'scale');
    assert(rgbScale.r >= 0 && rgbScale.r <= 255 && rgbScale.g >= 0 && rgbScale.g <= 255 && rgbScale.b >= 0 && rgbScale.b <= 255,
        'Масштабирование приводит все каналы в [0,255]');
    const rgbScale2 = Model.xyzToRgb(200, 50, 50, dataD65.xyz2rgb, 'scale');
    assert(rgbScale2.clipped === true,
        'Масштабирование помечает clipped=true при выходе за границы');
    return results;
}
class ColorApp {
    constructor() {
        this.testResultsDiv = document.getElementById('testResults');
        this.colorPreview = document.getElementById('colorPreview');
        this.illuminantSelect = document.getElementById('illuminantSelect');
        this.strategySelect = document.getElementById('clampStrategy');
        this.colorPicker = document.getElementById('colorPicker');
        this.state = {
            rgb: { r: 255, g: 0, b: 0 },
            xyz: { X: 41.24, Y: 21.26, Z: 1.93 },
            lab: { L: 53.24, a: 80.09, b: 67.20 },
            hsl: { h: 0, s: 1, l: 0.5 }
        };
        this.sliders = {};
        this.numberInputs = {};
        
        this.runTests();
        this.colorPreview.style.background = '#ff0000';
        
        document.querySelectorAll('.param').forEach(el => {
            const model = el.dataset.model;
            const comp = el.dataset.component;
            const slider = el.querySelector('input[type="range"]');
            const number = el.querySelector('input[type="number"]');
            const key = model + '.' + comp;
            this.sliders[key] = slider;
            this.numberInputs[key] = number;
            
            slider.addEventListener('input', () => {
                number.value = slider.value;
                this.onParamChange(model, comp, parseFloat(slider.value));
            });

            number.addEventListener('input', () => {
                let val = parseFloat(number.value);
                if (isNaN(val)) return;
                const min = parseFloat(slider.min);
                const max = parseFloat(slider.max);
                if (val < min) val = min;
                if (val > max) val = max;
                slider.value = val;
                this.onParamChange(model, comp, val);
            });
        });
        
        this.setFromRgb(255, 0, 0, 'init');
    }

    runTests() {
        const results = runTests();
        let output = '🔬 Результаты автотестов:\n';
        let allOk = true;
        results.forEach((res, i) => {
            const mark = res.ok ? 'OK' : 'BAD';
            output += `${mark} Тест ${i+1}: ${res.msg}\n`;
            if (!res.ok) allOk = false;
        });
        output += allOk ? '\nВсе тесты пройдены!' : '\nЕсть ошибки!';
        this.testResultsDiv.textContent = output;
    }

    fullUpdate() {
        const { r, g, b } = this.state.rgb;
        this.setFromRgb(r, g, b, 'fullUpdate');
    }

    setFromRgb(r, g, b, source, sourceClipped = false) {
        const illuminant = this.illuminantSelect.value;
        const strategy = this.strategySelect.value;
        const data = Model.getIlluminantData(illuminant);

        const xyz = Model.rgbToXyz(r, g, b, data.rgb2xyz);
        const lab = Model.xyzToLab(xyz.X, xyz.Y, xyz.Z, data.white);
        const rgbBack = Model.xyzToRgb(xyz.X, xyz.Y, xyz.Z, data.xyz2rgb, strategy);
        const hsl = Model.rgbToHsl(rgbBack.r, rgbBack.g, rgbBack.b);

        this.state.rgb = { r: rgbBack.r, g: rgbBack.g, b: rgbBack.b };
        this.state.xyz = xyz;
        this.state.lab = lab;
        this.state.hsl = hsl;
        this.lastClipped = sourceClipped || rgbBack.clipped;

        this.updateUI(source);
    }

    updateUI(source) {
        const { r, g, b } = this.state.rgb;
        const { X, Y, Z } = this.state.xyz;
        const { L, a, b: blab } = this.state.lab;
        const { h, s, l } = this.state.hsl;

        const setParam = (model, comp, val) => {
            const key = model + '.' + comp;
            if (this.sliders[key]) {
                this.sliders[key].value = val;
                this.numberInputs[key].value = val;
            }
        };
        setParam('xyz', 'X', X);
        setParam('xyz', 'Y', Y);
        setParam('xyz', 'Z', Z);
        setParam('lab', 'L', L);
        setParam('lab', 'a', a);
        setParam('lab', 'b', blab);
        setParam('hsl', 'H', h);
        setParam('hsl', 'S', s);
        setParam('hsl', 'L', l);

        if (source !== 'picker') {
            const hex = Model.rgbToHex(r, g, b);
            this.colorPicker.value = hex;
        }
        this.colorPreview.style.background = Model.rgbToHex(r, g, b);
    }

    onParamChange(model, comp, value) {
        const illuminant = this.illuminantSelect.value;
        const strategy = this.strategySelect.value;
        const data = Model.getIlluminantData(illuminant);

        let newRgb;
        if (model === 'xyz') {
            const x = (comp === 'X') ? value : this.state.xyz.X;
            const y = (comp === 'Y') ? value : this.state.xyz.Y;
            const z = (comp === 'Z') ? value : this.state.xyz.Z;
            const rgbRes = Model.xyzToRgb(x, y, z, data.xyz2rgb, strategy);
            newRgb = { r: rgbRes.r, g: rgbRes.g, b: rgbRes.b };
        } else if (model === 'lab') {
            const L = (comp === 'L') ? value : this.state.lab.L;
            const a = (comp === 'a') ? value : this.state.lab.a;
            const b = (comp === 'b') ? value : this.state.lab.b;
            const xyz2 = Model.labToXyz(L, a, b, data.white);
            const rgbRes = Model.xyzToRgb(xyz2.X, xyz2.Y, xyz2.Z, data.xyz2rgb, strategy);
            newRgb = { r: rgbRes.r, g: rgbRes.g, b: rgbRes.b };
        } else if (model === 'hsl') {
            const h = (comp === 'H') ? value : this.state.hsl.h;
            const s = (comp === 'S') ? value : this.state.hsl.s;
            const l = (comp === 'L') ? value : this.state.hsl.l;
            const rgb2 = Model.hslToRgb(h, s, l);
            newRgb = { r: rgb2.r, g: rgb2.g, b: rgb2.b };
        } else {
            return;
        }
        this.setFromRgb(newRgb.r, newRgb.g, newRgb.b, 'param', newRgb.clipped === true);
    }
}
document.addEventListener('DOMContentLoaded', () => {
    window.app = new ColorApp();
});