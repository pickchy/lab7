using System;
using System.Collections.Generic;
using System.Numerics;

class WaveletTransform
{
    private int SignalLength;
    private int DecompositionLevels;
    private double Pi = Math.PI;
    private List<Complex> LowPassFilter = new List<Complex>();
    private List<Complex> HighPassFilter = new List<Complex>();
    private List<List<Complex>> DetailBases = new List<List<Complex>>();
    private List<List<Complex>> ApproximationBases = new List<List<Complex>>();

    public WaveletTransform(int signalLength, int levels, int filterType)
    {
        SignalLength = signalLength;
        DecompositionLevels = levels;

 
        switch (filterType)
        {
            case 1: // Shannon 
                {
                    LowPassFilter.Add(1.0 / Math.Sqrt(2.0));
                    HighPassFilter.Add(1.0 / Math.Sqrt(2.0));
                    for (int i = 1; i < SignalLength; i++)
                    {
                        Complex temp = new Complex(
                            Math.Sqrt(2.0) / SignalLength * Math.Cos(Pi * i / SignalLength) * Math.Sin(Pi * i / 2.0) / Math.Sin(Pi * i / SignalLength),
                            -Math.Sqrt(2.0) / SignalLength * Math.Sin(Pi * i / SignalLength) * Math.Sin(Pi * i / 2.0) / Math.Sin(Pi * i / SignalLength)
                        );
                        LowPassFilter.Add(temp);
                        HighPassFilter.Add(Math.Pow(-1, i) * temp);
                    }
                }
                break;

            case 2: // Daubechies
                {
                    double a = 1.0 - Math.Sqrt(10.0);
                    double b = 1.0 + Math.Sqrt(10.0);
                    double c = Math.Sqrt(5.0 + 2.0 * Math.Sqrt(10.0));
                    double d = Math.Sqrt(2.0) / 32.0;

                    LowPassFilter.Add((b + c) * d);
                    LowPassFilter.Add((2.0 * a + 3.0 * b + 3.0 * c) * d);
                    LowPassFilter.Add((6.0 * a + 4.0 * b + 2.0 * c) * d);
                    LowPassFilter.Add((6.0 * a + 4.0 * b - 2.0 * c) * d);
                    LowPassFilter.Add((2.0 * a + 3.0 * b - 3.0 * c) * d);
                    LowPassFilter.Add((b - c) * d);

                    for (int i = 0; i < SignalLength - 6; i++)
                    {
                        LowPassFilter.Add(0);
                    }

                    HighPassFilter.Add(-LowPassFilter[1]);
                    HighPassFilter.Add(LowPassFilter[0]);
                    for (int i = 0; i < SignalLength - 6; i++)
                    {
                        HighPassFilter.Add(0);
                    }
                    HighPassFilter.Add(-LowPassFilter[5]);
                    HighPassFilter.Add(LowPassFilter[4]);
                    HighPassFilter.Add(-LowPassFilter[3]);
                    HighPassFilter.Add(LowPassFilter[2]);
                }
                break;

            default:
                throw new ArgumentException("Unknown filter type");
        }

        BuildBases();
    }

    private void BuildBases()
    {
        // Построение базисных функций
        List<List<Complex>> lowPassStages = new List<List<Complex>>() { new List<Complex>(LowPassFilter) };
        List<List<Complex>> highPassStages = new List<List<Complex>>() { new List<Complex>(HighPassFilter) };

        for (int level = 1; level < DecompositionLevels; level++)
        {
            int filterLength = SignalLength / (int)Math.Pow(2, level);
            lowPassStages.Add(new List<Complex>(new Complex[filterLength]));
            highPassStages.Add(new List<Complex>(new Complex[filterLength]));

            for (int j = 0; j < filterLength; j++)
            {
                int div = (int)Math.Pow(2, level);
                for (int k = 0; k < div; k++)
                {
                    lowPassStages[level][j] += lowPassStages[0][j + k * SignalLength / div];
                    highPassStages[level][j] += highPassStages[0][j + k * SignalLength / div];
                }
            }
        }

        DetailBases.Add(highPassStages[DecompositionLevels - 1]);
        ApproximationBases.Add(lowPassStages[DecompositionLevels - 1]);
    }

    public void Decompose(List<Complex> signal, int level, out List<Complex> approximationCoefficients, out List<Complex> detailCoefficients)
    {
        approximationCoefficients = new List<Complex>();
        detailCoefficients = new List<Complex>();

        int basisLength = LowPassFilter.Count / (int)Math.Pow(2, level);

        for (int i = 0; i < basisLength; i++)
        {
            detailCoefficients.Add(ComputeScalar(signal, DetailBases[i]));
            approximationCoefficients.Add(ComputeScalar(signal, ApproximationBases[i]));
        }
    }

    public void Reconstruct(List<Complex> approximationCoefficients, List<Complex> detailCoefficients, out List<Complex> reconstructedSignal)
    {
        reconstructedSignal = new List<Complex>(new Complex[SignalLength]);

        for (int i = 0; i < SignalLength; i++)
        {
            Complex approximation = 0;
            Complex detail = 0;

            for (int j = 0; j < approximationCoefficients.Count; j++)
            {
                approximation += approximationCoefficients[j] * ApproximationBases[j][i];
                detail += detailCoefficients[j] * DetailBases[j][i];
            }

            reconstructedSignal[i] = approximation + detail;
        }
    }

    private Complex ComputeScalar(List<Complex> vector1, List<Complex> vector2)
    {
        Complex result = 0;
        for (int i = 0; i < vector1.Count; i++)
        {
            result += vector1[i] * Complex.Conjugate(vector2[i]);
        }
        return result;
    }
}

class Program
{
    static double GenerateSignal(int n)
    {
        if ((100 <= n && n <= 150) || (200 <= n && n <= 250))
            return 1.0;
        if ((300 <= n && n <= 350) || (400 <= n && n <= 450))
            return 1.0 + 0.1 * Math.Cos(200.0 * Math.PI * n / 512.0);
        return 0.0;
    }

    static void Main()
    {
        List<Complex> signal = new List<Complex>();
        for (int i = 0; i < 512; i++)
        {
            signal.Add(new Complex(GenerateSignal(i), 0));
        }

        WaveletTransform wavelet = new WaveletTransform(512, 4, 2);

        wavelet.Decompose(signal, 4, out var approximation, out var detail);
        wavelet.Reconstruct(approximation, detail, out var reconstructedSignal);

        Console.WriteLine("Wavelet transform and signal reconstruction completed.");
    }
}
