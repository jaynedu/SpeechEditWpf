using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpeechEditWpf.signalprocess
{
    /// <summary>
    /// 复数类
    /// </summary>
    public class Complex
    {
        #region 属性
        private double real;      // 实部
        private double imag;      // 虚部
        public double Real { get { return real; } set { real = value; } }
        public double Imag { get { return imag; } set { imag = value; } }
        public double Mod { get { return Math.Sqrt(real * real + imag * imag); } }   // 模长
        public double Norm { get { return imag * imag + real * real; } }             // 范数
        #endregion
        #region 构造函数
        public Complex() { }

        public Complex(double real, double imag) 
        {
            this.real = real;
            this.imag = imag;
        }
        #endregion
        #region 运算符
        public static Complex operator +(Complex A, Complex B)
        {
            return new Complex(A.real + B.real, A.imag + B.imag);
        }
        public static Complex operator -(Complex A, Complex B)
        {
            return new Complex(A.real - B.real, A.imag - B.imag);
        }
        public static Complex operator *(Complex A, Complex B)
        {
            return new Complex(A.real * B.real - A.imag * B.imag, A.real * B.imag + A.imag * B.real);
        }
        public static Complex operator /(Complex A, Complex B)
        {
            if (A.real == 0 && B.imag == 0)
            {
                throw new Exception("除数不能为零！");
            }
            return new Complex((A * Conjugate(B)).real / B.Norm, (A * Conjugate(B)).imag / B.Norm);
        }
        #endregion
        #region 方法函数
        /// <summary>
        /// 共轭
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Complex Conjugate(Complex c)  
        {
            return new Complex(c.real, -c.imag);
        }

        /// <summary>
        /// 对复数序列取共轭
        /// </summary>
        /// <param name="c_sequence"></param>
        /// <returns></returns>
        public static Complex[] Conjugate(Complex[] c_sequence) 
        {
            Complex[] c_conj = new Complex[c_sequence.Length];
            for (int i = 0; i < c_sequence.Length; i++)
            {
                c_conj[i] = Conjugate(c_sequence[i]);
            }
            return c_conj;
        }

        /// <summary>
        /// 对复数矩阵取共轭
        /// </summary>
        /// <param name="c_matrix"></param>
        /// <returns></returns>
        public static Complex[][] Conjugate(Complex[][] c_matrix)
        {
            Complex[][] c_conj = new Complex[c_matrix.GetLength(0)][];
            for (int i = 0; i < c_matrix.GetLength(0); i++)
            {
                c_conj[i] = Conjugate(c_matrix[i]);
            }
            return c_conj;
        }

        /// <summary>
        /// 取复数序列每个元素的模
        /// </summary>
        /// <param name="c_sequence"></param>
        /// <returns></returns>
        public static double[] Abs(Complex[] c_sequence) 
        {
            double[] mod = new double[c_sequence.Length];
            for (int i = 0; i < c_sequence.Length; i++)
            {
                mod[i] = c_sequence[i].Mod;
            }
            return mod;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c_matrix"></param>
        /// <returns></returns>
        public static Matrix Abs(Complex[][] c_matrix)
        {
            double[][] mod = new double[c_matrix.GetLength(0)][];
            for (int i = 0; i < c_matrix.GetLength(0); i++)
            {
                mod[i] = Complex.Abs(c_matrix[i]);
            }
            return new Matrix(mod);
        }

        /// <summary>
        /// 取复数序列的实部
        /// </summary>
        /// <param name="c_sequence"></param>
        /// <returns></returns>
        public static double[] Reals(Complex[] c_sequence)
        {
            double[] reals = new double[c_sequence.Length];
            for (int i = 0; i < c_sequence.Length; i++)
            {
                reals[i] = c_sequence[i].real;
            }
            return reals;
        }

        /// <summary>
        /// 取复数矩阵的实部
        /// </summary>
        /// <param name="c_matrix"></param>
        /// <returns></returns>
        public static Matrix Reals(Complex[][] c_matrix) 
        {
            double[][] reals = new double[c_matrix.GetLength(0)][];
            for (int i = 0; i < c_matrix.GetLength(0); i++)
            {
                reals[i] = Complex.Reals(c_matrix[i]);
            }
            return new Matrix(reals);
        }

        /// <summary>
        /// 取复数序列的虚部
        /// </summary>
        /// <param name="c_sequence"></param>
        /// <returns></returns>
        public static double[] Imags(Complex[] c_sequence)
        {
            double[] imags = new double[c_sequence.Length];
            for (int i = 0; i < c_sequence.Length; i++)
            {
                imags[i] = c_sequence[i].imag;
            }
            return imags;
        }

        /// <summary>
        /// 倒谱分析
        /// </summary>
        /// <param name="c_matrix"></param>
        /// <returns></returns>
        public static Matrix Cepstrum(Complex[][] c_matrix) 
        {
            return Complex.Reals((Complex.Abs(c_matrix).Log(Math.E)).IFFT());
        }
        public override string ToString() 
        {
            if (this.real != 0 && this.imag != 0)
            {
                if (this.imag > 0)
                {
                    return string.Format("<{0}+j{1}>", this.real, this.imag);
                }
                else
                {
                    return string.Format("<{0}-j{1}>", this.real, -this.imag);
                }
            }
            else if (this.real == 0 && this.imag != 0)
            {
                if (this.imag > 0)
                {
                    return string.Format("<j{0}>", this.imag);
                }
                else
                {
                    return string.Format("<-j{1}>", this.real, -this.imag);
                }
            }
            else
            {
                return string.Format("<{0}>", this.real);
            }
        }
        #endregion
    }

    /// <summary>
    /// 序列类
    /// </summary>
    public class Vector : Operation
    {
        #region 属性
        private double[] values;      // 序列值
        public double[] Values { get { return values; } set { values = value; } }
        public int Length { get { return this.values.Length; } }  // 序列长度
        public double Sum { get { return this.values.Sum(); } } // 序列的和
        public double Min { get { return this.values.Min(); } } // 序列的最小值
        public double Max { get { return this.values.Max(); } } // 序列的最大值
        #endregion
        #region 构造函数
        public Vector(double[] vector)
        {
            this.values = vector;
        }

        public Vector() { }
        #endregion
        #region 运算符
        /// <summary>
        /// 序列相加
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Vector operator +(Vector A, Vector B)
        {
            if (A.Length != B.Length)
            {
                throw new Exception("输入参数维度不匹配");
            }
            else
            {
                double[] results = new double[A.Length];
                for (int i = 0; i < A.Length; i++)
                {
                    results[i] = A.values[i] + B.values[i];
                }
                return new Vector(results);
            }
        }

        /// <summary>
        /// 序列相减
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Vector operator -(Vector A, Vector B)
        {
            if (A.Length != B.Length)
            {
                throw new Exception("输入参数维度不匹配");
            }
            else
            {
                double[] results = new double[A.Length];
                for (int i = 0; i < A.Length; i++)
                {
                    results[i] = A.values[i] - B.values[i];
                }
                return new Vector(results);
            }
        }

        /// <summary>
        /// 序列相乘
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Vector operator *(Vector A, Vector B)
        {
            if (A.Length != B.Length)
            {
                throw new Exception("输入参数维度不匹配");
            }
            else
            {
                double[] results = new double[A.Length];
                for (int i = 0; i < A.Length; i++)
                {
                    results[i] = A.values[i] * B.values[i];
                }
                return new Vector(results);
            }
        }
        #endregion
        #region 方法函数
        public override string ToString()
        {
            string vec_str = "";
            foreach (double num in this.values)
            {
                vec_str += "<" + string.Format("{0:N2}", num) + "> ";
            }
            return vec_str;
        }

        /// <summary>
        /// 序列索引
        /// </summary>
        /// <param name="start"></param>
        /// <param name="stop"></param>
        /// <returns></returns>
        public Vector GetValues(int start, int stop)
        {
            if (start > stop)
            {
                throw new Exception("索引错误!");
            }
            else
            {
                int len = stop - start + 1;
                double[] v_slice = new double[len];
                for (int i = 0; i < len; i++)
                {
                    v_slice[i] = this.values[i + start];
                }
                return new Vector(v_slice);
            }
        }

        /// <summary>
        /// 序列后补零
        /// </summary>
        /// <param name="len_padded"></param>
        /// <param name="bits_num"></param>
        /// <returns></returns>
        public Vector ZeroPad(int len_padded, out int bits_num)
        {
            if (this.Length > len_padded)
            {
                throw new Exception("补零长度小于序列长度！");
            }
            bits_num = (int)Math.Ceiling(Math.Log(len_padded) / Math.Log(2));   // 计算补零后所需要的比特数
            double[] x = new double[len_padded];
            for (int i = 0; i < this.Length; i++)
            {
                x[i] = this.values[i];
            }

            return new Vector(x);
        }

        /// <summary>
        /// FFT运算
        /// </summary>
        /// <param name="N_fft"></param>
        /// <returns></returns>
        public Complex[] FFT(int N_fft)
        {
            /*
             *  <summary>
             *  功能: 快速傅里叶变换
             *  double[] x: 输入序列
             *  N_fft : FFT点数
             */

            double[] x_padded = this.ZeroPad(N_fft, out int bits_num).values;      // 后补零序列
            int[] inv_order = FFT_Order(N_fft, bits_num);                                    // 倒序下标
            Complex[] spectrum = new Complex[N_fft];                                     // 创建复数序列
            Complex[] butterfly_result = new Complex[2];                                 // 存放蝶形运算的中间结果
            int k, step, offset;

            // 将实数序列变为复数表示，并进行第一轮蝶形运算
            for (int i = 0; i < N_fft / 2; i++)
            {
                butterfly_result = Butterfly(new Complex(x_padded[inv_order[2 * i]], 0), new Complex(x_padded[inv_order[2 * i + 1]], 0), N_fft, 0);
                spectrum[2 * i] = butterfly_result[0];
                spectrum[2 * i + 1] = butterfly_result[1];
            }

            // 进行第fft_loop轮蝶形运算
            for (int fft_loop = bits_num - 2; fft_loop > 0; fft_loop--)
            {
                k = (int)Math.Pow(2, fft_loop);
                step = N_fft / 2 / k;
                offset = 2 * step;
                for (int i = 0; i < k; i++)
                {
                    for (int j = 0; j < step; j++)
                    {
                        butterfly_result = Butterfly(spectrum[i * offset + j], spectrum[i * offset + j + step], N_fft, j * k);
                        spectrum[i * offset + j] = butterfly_result[0];
                        spectrum[i * offset + j + step] = butterfly_result[1];
                    }
                }
            }

            // 进行最后一轮蝶形运算
            step = N_fft / 2;
            for (int j = 0; j < step; j++)
            {
                butterfly_result = Butterfly(spectrum[j], spectrum[j + step], N_fft, j);
                spectrum[j] = butterfly_result[0];
                spectrum[j + step] = butterfly_result[1];
            }
            return spectrum;
        }

        /// <summary>
        /// 语音分帧
        /// </summary>
        /// <param name="frame_len"></param>
        /// <param name="frame_overlap"></param>
        /// <returns></returns>
        public Matrix DivFrame(int frame_len, int frame_overlap)   // 
        {
            /*                                                    
            *  <summary>                                              
            *  功能：语音分帧                                             
            *  speech_sequence : 原始语音序列
            *  frame_num : 帧数
            */

            int step = frame_len - frame_overlap;
            int frame_num = (this.Length - frame_overlap) / step;               // 分出的帧数
            double[][] speech_divided = new double[frame_num][];             // 存储分帧后的数据矩阵

            for (int i = 0; i < frame_num; i++)       // 为每一行分配内存
            {
                speech_divided[i] = new double[frame_len];
            }

            for (int i = 0; i < frame_num - 1; i++)
            {
                for (int j = 0; j < frame_len; j++)
                {
                    speech_divided[i][j] = this.values[i * step + j];
                }
            }
            return new Matrix(speech_divided);
        }

        /// <summary>
        /// 对每个元素取对数
        /// </summary>
        /// <param name="e"></param>
        /// <returns></returns>
        public Vector Log(double e) 
        {
            double[] v_log = new double[this.Length];
            if (Math.E == e)
            {
                for (int i = 0; i < this.Length; i++)
                {
                    v_log[i] = Math.Log(this.values[i]);
                }
            }
            else if (e == 10)
            {
                for (int i = 0; i < this.Length; i++)
                {
                    v_log[i] = Math.Log10(this.values[i]);
                }
            }
            return new Vector(v_log);
        }

        /// <summary>
        /// 序列自相关分析
        /// </summary>
        /// <returns></returns>
        public Vector AutoCorrelation()
        {
            double[] Rn = new double[this.Length];
            for (int i = 0; i < this.Length; i++)
            {
                Rn[i] = (this.GetValues(0, this.Length - 1 - i) * this.GetValues(i, this.Length - 1)).Sum;
            }
            return new Vector(Rn);
        }
        #endregion
    }

    /// <summary>
    /// 矩阵类
    /// </summary>
    public class Matrix : Operation
    {
        #region 属性
        private double[][] values;      // 序列值
        public double[][] Values { get { return values; } set { values = value; } }
        public int Rows { get { return this.values.GetLength(0); } }  // 矩阵行数
        public int Cols { get { return this.values[0].Length; } }     //矩阵列数
        public double[] Min // 序列的最小值
        {
            get
            {
                double[] min = new double[this.Rows];
                for (int i = 0; i < this.Rows; i++)
                {
                    min[i] = this.values[i].Min();
                }
                return min;
            }
        }
        public double[] Max // 序列的最大值
        {
            get
            {
                double[] max = new double[this.Rows];
                for (int i = 0; i < this.Rows; i++)
                {
                    max[i] = this.values[i].Max();
                }
                return max;
            }
        }
        #endregion
        #region 构造函数
        public Matrix() { }

        public Matrix(double[][] mat)
        {
            this.values = mat;
        }
        #endregion
        #region 运算符
        /// <summary>
        /// 矩阵相加
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Matrix operator +(Matrix A, Matrix B)
        {
            if (A.Rows != B.Rows || A.Cols != B.Cols)
            {
                throw new Exception("输入参数维度不匹配");
            }
            else
            {
                double[][] results = new double[A.Rows][];
                for (int i = 0; i < A.Rows; i++)
                {
                    for (int j = 0; j < A.Cols; j++)
                    {
                        results[i][j] = A.values[i][j] + B.values[i][j];
                    }
                }
                return new Matrix(results);
            }
        }

        /// <summary>
        /// 矩阵相减
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Matrix operator -(Matrix A, Matrix B)
        {
            if (A.Rows != B.Rows || A.Cols != B.Cols)
            {
                throw new Exception("输入参数维度不匹配");
            }
            else
            {
                double[][] results = new double[A.Rows][];
                for (int i = 0; i < A.Rows; i++)
                {
                    for (int j = 0; j < A.Cols; j++)
                    {
                        results[i][j] = A.values[i][j] - B.values[i][j];
                    }
                }
                return new Matrix(results);
            }
        }

        /// <summary>
        /// 矩阵标量乘
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Matrix operator *(double A, Matrix B) 
        {
            for (int i = 0; i < B.Rows; i++)
            {
                for (int j = 0; j < B.Cols; j++)
                {
                    B.values[i][j] *= A;
                }
            }
            return B;
        }

        /// <summary>
        /// 矩阵乘
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Matrix operator *(Matrix A, Matrix B)
        {
            if (A.Cols != B.Rows)
            {
                throw new Exception("输入参数维度不匹配");
            }
            else
            {
                double[][] results = new double[A.Rows][];
                for (int i = 0; i < A.Rows; i++)
                {
                    results[i] = new double[B.Cols];
                    for (int j = 0; j < B.Cols; j++)
                    {
                        for (int k = 0; k < B.Rows; k++)
                        {
                            results[i][j] += A.values[i][k] * B.values[k][j];
                        }
                    }
                }
                return new Matrix(results);
            }
        }
        #endregion
        #region 方法函数
        public override string ToString()
        {
            string mat_str = "";
            for (int i = 0; i < this.Rows; i++)
            {
                for (int j = 0; j < this.Cols; j++)
                {
                    if (j == 0)
                    {
                        mat_str += "[" + string.Format("{0:N2}", this.values[i][j]) + " ";
                    }
                    else if (j == this.Cols - 1)
                    {
                        mat_str += string.Format("{0:N2}", this.values[i][j]) + "]\n";
                    }
                    else
                    {
                        mat_str += string.Format("{0:N2}", this.values[i][j]) + " ";
                    }
                }
            }
            return mat_str;
        }

        /// <summary>
        /// 矩阵转置
        /// </summary>
        /// <returns></returns>
        public Matrix Tranpose()
        {
            double[][] mat_tranposed = new double[this.Cols][];
            for (int i = 0; i < this.Cols; i++)
            {
                mat_tranposed[i] = new double[this.Rows];
                for (int j = 0; j < this.Rows; j++)
                {
                    mat_tranposed[i][j] = this.values[j][i];
                }
            }
            return new Matrix(mat_tranposed);
        }

        /// <summary>
        /// 对矩阵每个元素取底为e的对数
        /// </summary>
        /// <param name="e"></param>
        /// <returns></returns>
        public Matrix Log(double e)
        {
            double[][] mat_log = new double[this.Rows][];
            for (int i = 0; i < this.Rows; i++)
            {
                mat_log[i] = (new Vector(this.values[i]).Log(e)).Values;
            }
            return new Matrix(mat_log);
        }

        /// <summary>
        /// 矩阵每行后补零
        /// </summary>
        /// <param name="len_padded"></param>
        /// <param name="bits_num"></param>
        /// <returns></returns>
        public Matrix ZeroPad(int len_padded, out int bits_num)
        {
            if (this.Cols > len_padded)
            {
                throw new Exception("补零长度小于矩阵列数！");
            }

            bits_num = (int)Math.Ceiling(Math.Log(len_padded) / Math.Log(2));
            double[][] x = new double[this.Rows][];
            for (int i = 0; i < this.Rows; i++)
            {
                x[i] = new double[len_padded];
                for (int j = 0; j < this.Cols; j++)
                {
                    x[i][j] = this.values[i][j];
                }
            }
            return new Matrix(x);
        }

        /// <summary>
        /// 对矩阵每行进行FFT
        /// </summary>
        /// <param name="N_fft"></param>
        /// <returns></returns>
        public Complex[][] FFT(int N_fft)
        {
            double[][] x_padded = this.ZeroPad(N_fft, out int bits_num).values;     // 后补零序列
            int[] inv_order = FFT_Order(N_fft, bits_num);                                    // 倒序下标
            Complex[][] spectrum = new Complex[this.Rows][];                             // 创建复数序列
            Complex[] butterfly_result = new Complex[2];                                 // 存放蝶形运算的中间结果
            int k, step, offset;

            // 将实数序列变为复数表示，并进行第一轮蝶形运算
            for (int r = 0; r < this.Rows; r++)
            {
                spectrum[r] = new Complex[N_fft];
                for (int i = 0; i < N_fft / 2; i++)
                {
                    butterfly_result = Butterfly(new Complex(x_padded[r][inv_order[2 * i]], 0), new Complex(x_padded[r][inv_order[2 * i + 1]], 0), N_fft, 0);
                    spectrum[r][2 * i] = butterfly_result[0];
                    spectrum[r][2 * i + 1] = butterfly_result[1];
                }

                // 进行第fft_loop轮蝶形运算
                for (int fft_loop = bits_num - 2; fft_loop > 0; fft_loop--)
                {
                    k = (int)Math.Pow(2, fft_loop);
                    step = N_fft / 2 / k;
                    offset = 2 * step;
                    for (int i = 0; i < k; i++)
                    {
                        for (int j = 0; j < step; j++)
                        {
                            butterfly_result = Butterfly(spectrum[r][i * offset + j], spectrum[r][i * offset + j + step], N_fft, j * k);
                            spectrum[r][i * offset + j] = butterfly_result[0];
                            spectrum[r][i * offset + j + step] = butterfly_result[1];
                        }
                    }
                }

                // 进行最后一轮蝶形运算
                step = N_fft / 2;
                for (int j = 0; j < step; j++)
                {
                    butterfly_result = Butterfly(spectrum[r][j], spectrum[r][j + step], N_fft, j);
                    spectrum[r][j] = butterfly_result[0];
                    spectrum[r][j + step] = butterfly_result[1];
                }
            }
            return spectrum;
        }

        /// <summary>
        /// 对矩阵每行进行IFFT
        /// </summary>
        /// <returns></returns>
        public Complex[][] IFFT()
        {
            return Complex.Conjugate((((double)1 / this.Cols) * this).FFT(this.Cols));
        }

        /// <summary>
        /// 短时能量分析
        /// </summary>
        /// <returns></returns>
        public double[] STEnergy()
        {
            double[] En = new double[this.Rows];
            Matrix mat = this * this.Tranpose();
            for (int i = 0; i < this.Rows; i++)
            {
                En[i] += mat.values[i][i];
            }
            return En;
        }

        /// <summary>
        /// 短时过零率分析
        /// </summary>
        /// <returns></returns>
        public double[] ZCR()         // 
        {
            double[] Zn = new double[this.Rows];
            for (int i = 0; i < this.Rows; i++)
            {
                for (int j = 1; j < this.Cols; j++)
                {
                    Zn[i] += Math.Abs(Math.Sign(this.values[i][j]) - Math.Sign(this.values[i][j - 1])) / 2;
                }
            }
            return Zn;
        }

        /// <summary>
        /// 矩阵每行数据进行自相关分析
        /// </summary>
        /// <returns></returns>
        public Matrix AutoCorrelation() 
        {
            double[][] Rn = new double[this.Rows][];
            for (int i = 0; i < this.Rows; i++)
            {
                Rn[i] = (new Vector(this.Values[i]).AutoCorrelation()).Values;
            }
            return new Matrix(Rn);
        }
        #endregion
    }

    /// <summary>
    /// WAV文件结构内容描述
    /// </summary>
    public struct WavInfo
    {
        public int SamplePerSec;    // 采样率
        public string Channels;     // 通道数信息
        public string BitsPerSample;
        public int Data_Size;
        public int Data_Len;
        public double Wav_Time;     // wav数据时间长度
    }
}
