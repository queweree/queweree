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