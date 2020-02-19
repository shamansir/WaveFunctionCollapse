/*
The MIT License(MIT)
Copyright(c) mxgmn 2016.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software.
*/

using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Collections.Generic;

class OverlappingModel : Model
{
    int N; // `N` is the number of pixels in the side of the pattern
    byte[][] patterns; // an array of unique patterns found in the source image
    List<Color> colors; // an array of all possible colors found in the source image
    int ground; // the ground value, should not be higher than image height

    public OverlappingModel(string name, int N, int width, int height, bool periodicInput, bool periodicOutput, int symmetry, int ground)
        : base(width, height)
    {
        this.N = N; // store the `N` (the number of pixels in the side of the pattern)
        // `periodicOutput` means that we treat the output image as infinitely
        // repeating in every direction and so there are no borderline
        // pixels: every edge wraps and gets connected to the opposite edge;
        // don't mix with the `periodicInput` which means the same but
        // for the source image, and so is called `periodicInput` versus
        // `periodic` below;
        periodic = periodicOutput;

        // SMX and SMY are pixel dimensions of the source image
        var bitmap = new Bitmap($"samples/{name}.png");
        int SMX = bitmap.Width, SMY = bitmap.Height;
        // this `SMXxSMY` will store the indices of unique colors ,
        // _not_ RGBA or the color value in any sense, but the color ID
        // for every pixel
        byte[,] sample = new byte[SMX, SMY];
        colors = new List<Color>();

        // collect all the unique colors from the image,
        // store them into the `colors` list;
        for (int y = 0; y < SMY; y++) for (int x = 0; x < SMX; x++)
            {
                Color color = bitmap.GetPixel(x, y);

                int i = 0;
                foreach (var c in colors)
                {
                    if (c == color) break;
                    i++;
                }

                if (i == colors.Count) colors.Add(color);
                sample[x, y] = (byte)i;
            }

        // `C` is now the number of unique colors in the image
        int C = colors.Count;
        // and `W` is `C ^ (N ^ 2)`, so the number of
        // all the possible x <-> y <-> (unique color index) combinations:
        // for every pixel for any pattern of size NxN there is
        // an equal possibility for one of the unique colors
        // to be present at that place;
        long W = C.ToPower(N * N);

        // build a pattern:
        // take some function which returns color value
        // for a given position (`x` & `y`), and build a
        // flat array of the size (`N * N`) filled with
        // the values returned from such function;
        byte[] pattern(Func<int, int, byte> f)
        {
            byte[] result = new byte[N * N];
            for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) result[x + y * N] = f(x, y);
            return result;
        };

        // so this function gets the pattern from the source by its location
        // (`x` & `y`) and returns the flat array of its pixels;
        byte[] patternFromSample(int x, int y) => pattern((dx, dy) => sample[(x + dx) % SMX, (y + dy) % SMY]);
        // and this one rotates the given pattern by 90 degrees and returns the rearranged pixels
        byte[] rotate(byte[] p) => pattern((x, y) => p[N - 1 - y + x * N]);
        // and this one mirrors the given pattern and returns the rearranged pixels
        byte[] reflect(byte[] p) => pattern((x, y) => p[N - 1 - x + y * N]);

        // `index` here is like the unique hash sum (fingerprint) of the pattern,
        // so for the patterns with the same color IDs at the same positions,
        // their hash sums (here, indices) are equal
        long index(byte[] p)
        {
            long result = 0, power = 1;
            for (int i = 0; i < p.Length; i++)
            {
                // the `i` used in the reverse order because
                // the `unpacking` process uses the subtraction
                // to restore the color IDs
                result += p[p.Length - 1 - i] * power;
                // remember that C is also the amount
                // of all the possibilities for the specific pixel
                // and the number of unique colors in the source image
                power *= C;
            }
            return result;
        };

        // get all the pixels of the pattern of size `N` using its known
        // hash sum, which is revertible
        byte[] patternFromIndex(long ind)
        {
            long residue = ind, power = W;
            byte[] result = new byte[N * N];

            for (int i = 0; i < result.Length; i++)
            {
                power /= C;
                int count = 0;

                while (residue >= power)
                {
                    residue -= power;
                    count++;
                }

                result[i] = (byte)count;
            }

            return result;
        };

        // the dictionary, mapping the unique hash sums of the patterns
        // to their weights, where `weight` is the amount of times the pattern
        // occured in the source, considering its reflections and rotations
        Dictionary<long, int> weights = new Dictionary<long, int>();
        // the list of all the hash sums of the unique patterns
        List<long> ordering = new List<long>();

        // collect the unique patterns in the source image,
        // by visiting every pixel in that image, extracting
        // the pattern of size `NxN` starting from that pixel,
        // rotating and reflecting it and then comparing the hash
        // sum of every such modification to what was found before
        // — and if that hash sum already was met, increasing the `weight`
        // value for that hash sum;
        // if `periodicInput` is disabled, stop before the patterns
        // would overflow image borders;
        for (int y = 0; y < (periodicInput ? SMY : SMY - N + 1); y++) for (int x = 0; x < (periodicInput ? SMX : SMX - N + 1); x++)
            {
                // the array to store all the symmetries of the pattern
                byte[][] ps = new byte[8][];

                // extract the pattern at that position and rotate and reflect it
                // in four directions
                // (so four directions + four reflections == eight symmetries)
                ps[0] = patternFromSample(x, y);
                ps[1] = reflect(ps[0]);
                ps[2] = rotate(ps[0]);
                ps[3] = reflect(ps[2]);
                ps[4] = rotate(ps[2]);
                ps[5] = reflect(ps[4]);
                ps[6] = rotate(ps[4]);
                ps[7] = reflect(ps[6]);

                for (int k = 0; k < symmetry; k++)
                {
                    // calculate the hash sum for that symmetry
                    long ind = index(ps[k]);
                    // increase the weight for that hash sum
                    // if it was already met
                    if (weights.ContainsKey(ind)) weights[ind]++;
                    else
                    {
                        // or instead add it to the storage with the weight of `1`
                        weights.Add(ind, 1);
                        ordering.Add(ind);
                    }
                }
            }

        // `T` is the number of unique patterns found
        T = weights.Count;
        // make `ground` to not be higher than the number of unique patterns
        this.ground = (ground + T) % T;
        // prepare the `patterns` array to have the size of `T`
        patterns = new byte[T][];
        base.weights = new double[T];

        int counter = 0;
        // walk through all the unique patterns hash-sums (in the order they were
        // discovered) and restore the patterns back to the matrices of the color IDs,
        // and also fill in the `base` class storage of weights with the `weights`
        // calculated before
        foreach (long w in ordering)
        {
            patterns[counter] = patternFromIndex(w);
            base.weights[counter] = weights[w];
            counter++;
        }

        // given the `p1` as the pattern of size `N*N` and
        // `p2` as the source image of whichever size,
        // take the `NxN` rectangle in the source image at the
        // position `(dx, dy)` and compare all the values
        // with the pattern, return `true` when they are totally equal
        bool agrees(byte[] p1, byte[] p2, int dx, int dy)
        {
            int xmin = dx < 0 ? 0 : dx, xmax = dx < 0 ? dx + N : N, ymin = dy < 0 ? 0 : dy, ymax = dy < 0 ? dy + N : N;
            for (int y = ymin; y < ymax; y++) for (int x = xmin; x < xmax; x++) if (p1[x + N * y] != p2[x - dx + N * (y - dy)]) return false;
            return true;
        };

        // `DX` = { -1, 0, 1, 0 }
        // `DY` = { 0, 1, 0, -1 }
        // for the four cases, for every unique pattern,
        // propagator stores the indices of the unique patterns that
        // match with that pattern at the position of `(DX[case_index], DY[case_index])`
        propagator = new int[4][][];
        for (int d = 0; d < 4; d++)
        {
            // `T` is the number of unique patterns
            propagator[d] = new int[T][];
            for (int t = 0; t < T; t++)
            {
                List<int> list = new List<int>();
                // for every unique pattern, take the other unique pattern
                // from the same list and if they agree (match) at the
                // `(DX[case_index], DY[case_index])` point,
                // add the latter to the list of successes at this position
                for (int t2 = 0; t2 < T; t2++) if (agrees(patterns[t], patterns[t2], DX[d], DY[d])) list.Add(t2);
                propagator[d][t] = new int[list.Count];
                for (int c = 0; c < list.Count; c++) propagator[d][t][c] = list[c];
            }
        }

        // `propagator` actually stores the information about how the patterns match each other,
        // (independently of the source) at different offsets:
        // `(-1, 0)`, `(0, 1)`, `(1, 0)` and `(0, -1)`.
        // notice that it's only in four directions: `N`, `W`, `S`, `E`
    }

    protected override bool OnBoundary(int x, int y) => !periodic && (x + N > FMX || y + N > FMY || x < 0 || y < 0);

    public override Bitmap Graphics()
    {
        Bitmap result = new Bitmap(FMX, FMY);
        int[] bitmapData = new int[result.Height * result.Width];

        if (observed != null)
        {
            for (int y = 0; y < FMY; y++)
            {
                int dy = y < FMY - N + 1 ? 0 : N - 1;
                for (int x = 0; x < FMX; x++)
                {
                    int dx = x < FMX - N + 1 ? 0 : N - 1;
                    Color c = colors[patterns[observed[x - dx + (y - dy) * FMX]][dx + dy * N]];
                    bitmapData[x + y * FMX] = unchecked((int)0xff000000 | (c.R << 16) | (c.G << 8) | c.B);
                }
            }
        }
        else
        {
            for (int i = 0; i < wave.Length; i++)
            {
                int contributors = 0, r = 0, g = 0, b = 0;
                int x = i % FMX, y = i / FMX;

                for (int dy = 0; dy < N; dy++) for (int dx = 0; dx < N; dx++)
                    {
                        int sx = x - dx;
                        if (sx < 0) sx += FMX;

                        int sy = y - dy;
                        if (sy < 0) sy += FMY;

                        int s = sx + sy * FMX;
                        if (OnBoundary(sx, sy)) continue;
                        for (int t = 0; t < T; t++) if (wave[s][t])
                            {
                                contributors++;
                                Color color = colors[patterns[t][dx + dy * N]];
                                r += color.R;
                                g += color.G;
                                b += color.B;
                            }
                    }

                bitmapData[i] = unchecked((int)0xff000000 | ((r / contributors) << 16) | ((g / contributors) << 8) | b / contributors);
            }
        }

        var bits = result.LockBits(new Rectangle(0, 0, result.Width, result.Height), ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);
        System.Runtime.InteropServices.Marshal.Copy(bitmapData, 0, bits.Scan0, bitmapData.Length);
        result.UnlockBits(bits);

        return result;
    }

    protected override void Clear()
    {
        base.Clear();

        if (ground != 0)
        {
            for (int x = 0; x < FMX; x++)
            {
                for (int t = 0; t < T; t++) if (t != ground) Ban(x + (FMY - 1) * FMX, t);
                for (int y = 0; y < FMY - 1; y++) Ban(x + y * FMX, ground);
            }

            Propagate();
        }
    }
}
