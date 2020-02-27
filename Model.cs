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
        // those four numbers represent the amount of tiles/patterns compatible
        // with this tile/pattern for every direction;
        // this is the state of propagation, described in details below
        compatible = new int[wave.Length][][];

        // preparing the `wave` and `compatible` arrays
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

        // `sumOfOnes` stores the amount of possibilities left to consider,
        // for every slot in the output
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
        // the most possible minimum value
        double min = 1E+3;
        int argmin = -1;

        // the main calculation of the entropy,
        // `wave.Length` is equal to `FMX * FMY` here
        for (int i = 0; i < wave.Length; i++)
        {
            // `OnBoundary` is overloaded in the children classes
            // and skips (doesn't affect) the cells which are close to
            // the edges; for `OverlappingModel` it is only if output
            // is not `periodic` and if position falls into the pattern on
            // the edge; for `SimpleTiledModel` it just controls if
            // position is not outside the edges, also only if output is not
            // `periodic`; `i % FMX` and `i / FMX` are `x` and `y`;
            if (OnBoundary(i % FMX, i / FMX)) continue;

            // as it was mentioned above, `sumsOfOnes` holds the amount
            // of possibilities left to consider for the output slot at the specified position
            int amount = sumsOfOnes[i];
            // if it is turned out that there are no more options at any position,
            // then it's a contradictions and the algorithm should completely stop;
            if (amount == 0) return false;

            double entropy = entropies[i];
            // if there are still some options to consider (if amount equals to 1,
            // the cell is solved) and entropy is less than anything we've checked
            // before...
            if (amount > 1 && entropy <= min)
            {
                // calculate the very small random noise value, to avoid sticking
                // to the same solution over time
                double noise = 1E-6 * random.NextDouble();
                // and if entropy, with some random noise, is still minimal...
                if (entropy + noise < min)
                {
                    // update the minimum entropy value to this entropy (with the noise)
                    min = entropy + noise;
                    // set `argmin` flag to `i`, meaning that we did found at least one minimum
                    // entropy value and the latest one was located at position `i`
                    argmin = i;
                }
            }
        }

        // if no minimum entropy was found, that means that the state is observed:
        // all the solutions are known
        if (argmin == -1)
        {
            observed = new int[FMX * FMY];
            // fill in the state of observation with the indices of the tiles/patterns for
            // which `wave` had `true` in that position, one index per position,
            // and it is the first such index found in the list of possibilities (notice `break`)
            for (int i = 0; i < wave.Length; i++) for (int t = 0; t < T; t++) if (wave[i][t]) { observed[i] = t; break; }
            // now that `observed` array contains all the solutions, exit the observation
            return true;
        }

        // `argmin` holds the index of the output slot with the minimim entropy;

        // fill in the `distribution` array (it has the size equal to the number of unique
        // tiles/patterns found in the input) with the weights values for every tile/pattern
        // for which `wave` has `true` value or `0` otherwise
        double[] distribution = new double[T];
        for (int t = 0; t < T; t++) distribution[t] = wave[argmin][t] ? weights[t] : 0;
        // take index of the random tile from the distribution array and store it in `r` variable
        int r = distribution.Random(random.NextDouble());
        // ?? seems the values in the `distribution` array are not used in any matter

        // from `wave`, get the list of boolean values for the possibilities left
        // for the cell with the minimum entropy
        bool[] w = wave[argmin];
        // for all the possibilities, if the value in the `wave` for that possibility is
        // falsey, or the tile/pattern index by accident is equal to the random `r` index,
        // ban that possibility (set the value in the `wave` to `false`)
        for (int t = 0; t < T; t++) if (w[t] != (t == r)) Ban(argmin, t);

        return null;
    }

    protected void Propagate()
    {
        // stack contains the banned possibilities (tiles/patterns)
        // at particular positions in the output;
        // the stack can grow from inside this loop, so
        while (stacksize > 0)
        {
            // pop the next pair of position and possibility from `stack`
            var e1 = stack[stacksize - 1];
            stacksize--;

            int i1 = e1.Item1;
            // convert index to the X/Y coordinates back
            int x1 = i1 % FMX, y1 = i1 / FMX;

            // for every direction
            for (int d = 0; d < 4; d++)
            {
                // find `x2` and `y2` which are `x1` and `y1` moved
                // in the given direction
                int dx = DX[d], dy = DY[d];
                int x2 = x1 + dx, y2 = y1 + dy;
                // if the new position turned out to be on boundary, skip it
                if (OnBoundary(x2, y2)) continue;

                // wrap the edges of the output, so the left edge and right edge
                // are bound together and top edge and bottom edges as well
                if (x2 < 0) x2 += FMX;
                else if (x2 >= FMX) x2 -= FMX;
                if (y2 < 0) y2 += FMY;
                else if (y2 >= FMY) y2 -= FMY;

                // calculate the index from this position
                int i2 = x2 + y2 * FMX;
                // `e1.Item2` is the index of the of the tile/pattern which was banned
                // during the latest observation;
                // `propagator` contains the information about the tiles/patterns that
                // match to the given tile/pattern (the banned one, in this case)
                // at every direction;
                // so here we load the list of the tiles/patterns that match the banned pattern
                // at the direction currently under consideration;
                int[] p = propagator[d][e1.Item2];
                // `compatible` array stores the number of the tiles/patterns matching to that
                // tile/pattern in the output, but unlike `sumOfOnes`, this one is about the
                // neigbourgs of the slot rather that the slot itself
                int[][] compat = compatible[i2];

                // for every tile/pattern that were matching the banned tile at this position
                for (int l = 0; l < p.Length; l++)
                {
                    // store its ID in `t2`
                    int t2 = p[l];

                    // by this ID, find how much there were compatible tiles/patterns left for
                    // this tile/pattern
                    int[] comp = compat[t2];

                    // reduce the number of such tiles/patterns
                    comp[d]--;
                    // and if it became zero, ban this possibility as well;
                    // NB: that action puts the newly banned possibility in stack
                    // and so increases the `stacksize` again, so provocating (or propagating)
                    // that action again and again until all the questionable combinations
                    // are solved;
                    if (comp[d] == 0) Ban(i2, t2);
                }
            }
        }
    }

    public bool Run(int seed, int limit)
    {
        if (wave == null) Init();

        // clear is not actually clearing values, but sets the
        // values to initial; including setting all `wave` values to `true`
        // and loading values from `propagator` (overloaded in the child model)
        // into the `compatible` array
        Clear();
        random = new Random(seed);

        // try until the limit of tries is reached
        for (int l = 0; l < limit || limit == 0; l++)
        {
            // observe: if result of the observation is `null` (observations are done
            // without explicit result), do the propogation and continue; if result is
            // `true`, then the output is fully observed, so exit with a success;
            // if result is `false` — that means the output cannot be solved,
            // it's a contracdiction, exit with failure;
            bool? result = Observe();
            if (result != null) return (bool)result;
            Propagate();
        }

        return true;
    }

    protected void Ban(int i, int t)
    {
        // banning is setting the exact possibility in the given position to `false`;
        // `i` is the index of the X/Y position, but in the flatten array, `t` is the
        // index of the tile/pattern which would fit at this position
        wave[i][t] = false;

        // clear the `compatible` counter for that position and tile/pattern to `0`
        // in all four directions
        int[] comp = compatible[i][t];
        for (int d = 0; d < 4; d++) comp[d] = 0;
        // we put the pair of the X/Y position and the posiibility (tile/pattern) index
        // onto the stack
        stack[stacksize] = (i, t);
        stacksize++;

        // reduce the amount of tiles/patterns which can fit at this position
        sumsOfOnes[i] -= 1;
        // reduce the sum of weights by the weight of this tile/pattern for this position
        // in the output
        sumsOfWeights[i] -= weights[t];
        // as well as logarithm of weight
        sumsOfWeightLogWeights[i] -= weightLogWeights[t];

        // recalculate the entropy for this position
        double sum = sumsOfWeights[i];
        entropies[i] = Math.Log(sum) - sumsOfWeightLogWeights[i] / sum;
    }

    protected virtual void Clear()
    {
        // through all the output positions
        for (int i = 0; i < wave.Length; i++)
        {
            // and every unique tile/pattern
            for (int t = 0; t < T; t++)
            {
                // enable the possibility of it to appeear to `true` in the `wave` array
                wave[i][t] = true;
                // collect the number of `compatible` tiles/patterns for that position
                // when it is shifted in given direction, by using the data from the `propagator`
                // but in opposite direction;
                for (int d = 0; d < 4; d++) compatible[i][t][d] = propagator[opposite[d]][t].Length;
            }

            // `sumOfOnes` contains the amount of possible answers, which is equal to
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
