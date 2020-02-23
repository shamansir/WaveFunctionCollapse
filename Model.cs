/*
The MIT License(MIT)
Copyright(c) mxgmn 2016.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software.
*/

using System;

abstract class Model
{
    protected bool[][] wave;

    // `propagator` stores the information
    // about how each tile/pattern matches each other pattern
    // (or itself) in four directions: `N`, `W`, `S`, `E`;
    // so, for every direction, for every tile, there's a list
    // of matching tiles for this particular tile in
    // that particular direction
    protected int[][][] propagator;
    int[][][] compatible;
    protected int[] observed;

    (int, int)[] stack;
    int stacksize;

    protected Random random;
    // `FMXxFMY` is the size of the output image
    // `T` is the number of unique tiles or patterns
    protected int FMX, FMY, T;
    protected bool periodic;

    protected double[] weights;
    double[] weightLogWeights;

    int[] sumsOfOnes;
    double sumOfWeights, sumOfWeightLogWeights, startingEntropy;
    double[] sumsOfWeights, sumsOfWeightLogWeights, entropies;

    protected Model(int width, int height)
    {
        FMX = width;
        FMY = height;
    }

    void Init()
    {
        // initialize the `wave` as the flat array where every item is representing
        // the "slot" in the ouput image, where _slot_ for `OverlappingModel` represents the
        // single pixel, and for `SimpleTiled` model is the tile to place;
        // the value of every such item is another array of the size `T`
        // (number of unique "tiles"), so it contains the boolean value for every possible outcome: // if the outcoume still possible for that particular slot, or already not;
        wave = new bool[FMX * FMY][];
        // this array is almost of the same shape (every _slot_ of the ouput image bound
        // to the list of the possible outcomes), but contains arrays of four integers
        // instead of one boolean, so each outcome has a corresponding four-integer array
        // associated with it, four is the number of possible directions: `N`, `S`, `W`, `E`;
        // this is the state of propagation, described in details below
        compatible = new int[wave.Length][][];

        // preparing and filling up the `wave` and `compatible` arrays
        for (int i = 0; i < wave.Length; i++)
        {
            wave[i] = new bool[T];
            compatible[i] = new int[T][];
            for (int t = 0; t < T; t++) compatible[i][t] = new int[4];
        }

        weightLogWeights = new double[T];
        sumOfWeights = 0;
        sumOfWeightLogWeights = 0;

        // use the values of `weights` which were calculated in the `Init` method of
        // a child class (`OverlappingModel` or `SimpleTiledModel`): `T` is the number of
        // unique tiles/patterns found in the source, and the value of `weights[t]`
        // is the number of times this tile/pattern appeared in the source
        for (int t = 0; t < T; t++)
        {
            weightLogWeights[t] = weights[t] * Math.Log(weights[t]);
            sumOfWeights += weights[t];
            sumOfWeightLogWeights += weightLogWeights[t];
        }

        startingEntropy = Math.Log(sumOfWeights) - sumOfWeightLogWeights / sumOfWeights;

        // for the sake of caching and to not repeat unnecessary calculations,
        // all the precalculated values are stored in the corresponding arrays
        // of the size equal to the output image size, being flatten
        sumsOfOnes = new int[FMX * FMY];
        sumsOfWeights = new double[FMX * FMY];
        sumsOfWeightLogWeights = new double[FMX * FMY];
        entropies = new double[FMX * FMY];

        // `stack`, where each element is the coordinate (TODO),
        // and its size is the output image size multiplied
        // by the number of unique patterns/tiles
        stack = new (int, int)[wave.Length * T];
        stacksize = 0;
    }

    bool? Observe()
    {
        double min = 1E+3;
        int argmin = -1;

        for (int i = 0; i < wave.Length; i++)
        {
            if (OnBoundary(i % FMX, i / FMX)) continue;

            int amount = sumsOfOnes[i];
            if (amount == 0) return false;

            double entropy = entropies[i];
            if (amount > 1 && entropy <= min)
            {
                double noise = 1E-6 * random.NextDouble();
                if (entropy + noise < min)
                {
                    min = entropy + noise;
                    argmin = i;
                }
            }
        }

        if (argmin == -1)
        {
            observed = new int[FMX * FMY];
            for (int i = 0; i < wave.Length; i++) for (int t = 0; t < T; t++) if (wave[i][t]) { observed[i] = t; break; }
            return true;
        }

        double[] distribution = new double[T];
        for (int t = 0; t < T; t++) distribution[t] = wave[argmin][t] ? weights[t] : 0;
        int r = distribution.Random(random.NextDouble());

        bool[] w = wave[argmin];
        for (int t = 0; t < T; t++) if (w[t] != (t == r)) Ban(argmin, t);

        return null;
    }

    protected void Propagate()
    {
        while (stacksize > 0)
        {
            var e1 = stack[stacksize - 1];
            stacksize--;

            int i1 = e1.Item1;
            int x1 = i1 % FMX, y1 = i1 / FMX;

            for (int d = 0; d < 4; d++)
            {
                int dx = DX[d], dy = DY[d];
                int x2 = x1 + dx, y2 = y1 + dy;
                if (OnBoundary(x2, y2)) continue;

                if (x2 < 0) x2 += FMX;
                else if (x2 >= FMX) x2 -= FMX;
                if (y2 < 0) y2 += FMY;
                else if (y2 >= FMY) y2 -= FMY;

                int i2 = x2 + y2 * FMX;
                int[] p = propagator[d][e1.Item2];
                int[][] compat = compatible[i2];

                for (int l = 0; l < p.Length; l++)
                {
                    int t2 = p[l];
                    int[] comp = compat[t2];

                    comp[d]--;
                    if (comp[d] == 0) Ban(i2, t2);
                }
            }
        }
    }

    public bool Run(int seed, int limit)
    {
        if (wave == null) Init();

        Clear();
        random = new Random(seed);

        for (int l = 0; l < limit || limit == 0; l++)
        {
            bool? result = Observe();
            if (result != null) return (bool)result;
            Propagate();
        }

        return true;
    }

    protected void Ban(int i, int t)
    {
        wave[i][t] = false;

        int[] comp = compatible[i][t];
        for (int d = 0; d < 4; d++) comp[d] = 0;
        stack[stacksize] = (i, t);
        stacksize++;

        sumsOfOnes[i] -= 1;
        sumsOfWeights[i] -= weights[t];
        sumsOfWeightLogWeights[i] -= weightLogWeights[t];

        double sum = sumsOfWeights[i];
        entropies[i] = Math.Log(sum) - sumsOfWeightLogWeights[i] / sum;
    }

    protected virtual void Clear()
    {
        for (int i = 0; i < wave.Length; i++)
        {
            for (int t = 0; t < T; t++)
            {
                wave[i][t] = true;
                for (int d = 0; d < 4; d++) compatible[i][t][d] = propagator[opposite[d]][t].Length;
            }

            sumsOfOnes[i] = weights.Length;
            sumsOfWeights[i] = sumOfWeights;
            sumsOfWeightLogWeights[i] = sumOfWeightLogWeights;
            entropies[i] = startingEntropy;
        }
    }

    protected abstract bool OnBoundary(int x, int y);
    public abstract System.Drawing.Bitmap Graphics();

    protected static int[] DX = { -1, 0, 1, 0 };
    protected static int[] DY = { 0, 1, 0, -1 };
    static int[] opposite = { 2, 3, 0, 1 };
}
