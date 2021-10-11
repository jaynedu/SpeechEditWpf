using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpeechEditWpf.signalprocess
{
    public class Operation
    {
        /// <summary>
        /// 蝶形运算
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="N_fft">做FFT的点数</param>
        /// <param name="k">第 k 轮蝶形运算</param>
        /// <returns></returns>
        public Complex[] Butterfly(Complex X, Complex Y, int N_fft, int k)
        {
            Complex Wnk = new Complex(Math.Cos(2 * Math.PI * k / (double)N_fft), -Math.Sin(2 * Math.PI * k / (double)N_fft));     // 旋转因子
            Complex[] Output = new Complex[2];       // 蝶形运算的输出                                                                                               
            Output[0] = X + Y * Wnk;
            Output[1] = X - Y * Wnk;
            return Output;
        }

        /// <summary>
        /// 根据输入数据的的序号的二进制倒转来重新排序
        /// </summary>
        /// <param name="N_fft"></param>
        /// <param name="bits_num"></param>
        /// <returns>返回值y为行向量</returns>
        /// 输入x(0) x(1) x(2) x(3) x(4) x(5) x(6) x(7)
        /// 二进制序号000 001 010 011 100 101 110 111
        /// 二进制反序000 100 010 110 001 101 011 111
        /// 输出x(0) x(4) x(2) x(6) x(1) x(5) x(3) x(7)
        public int[] FFT_Order(int N_fft, int bits_num) 
        {

            int[] inv_index = new int[N_fft];
            char[] bisequence;
            for (int i = 0; i < N_fft; i++)
            {
                bisequence = Convert.ToString(i, 2).PadLeft(bits_num, '0').ToArray();
                Array.Reverse(bisequence, 0, bits_num);
                inv_index[i] = Convert.ToInt32(new string(bisequence), 2);
            }
            return inv_index;
        }
    }
}
