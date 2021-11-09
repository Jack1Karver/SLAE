using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SLAE
{
    class Program
    {
        static double[,] matrixConst = { { 3.43, 3.38, 3.09 }, { 4.17, 4, 3.65 }, { 4.30, 4.10, 3.67 } };
                
        static double[] fmembConst = { 5.52, 6.93, 7.29 };
        static int size = fmembConst.Length;
        static double[] fmemb = new double[size]; 
        static double[] y = new double[size];
        static double[,] matrix = new double[size, size];
        static double[,] L = new double[size, size];
        static double[,] L1 = new double[size,size];        
        static double[] res = new double[size];
        static void Main(string[] args)
        {
            while (true) {
                Console.WriteLine("\nSelect method:\n" +
                    "1 -- Gauss Method\n" + "2 -- Cholesky method");
                switch (Console.ReadLine()) {
                    case "1":
                        {
                            Gauss();
                            break;
                        }
                    case "2":
                        {
                            Cholesky();
                            break;
                        }

            }
               
            }
            
        }

        public static void Gauss() {           

            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    L[i, j] = matrixConst[i, j];                    
                }
                fmemb[i] = fmembConst[i];               
            }
            WriteMatrix(L, fmemb);
            for (int i = 0; i < size-1; i++)
            {
                ReFactor(i);
                
                Division(i);                
                Substraction(i);                
            }
            fmemb[size - 1] /= L[size - 1, size - 1];
            L[size - 1, size - 1] /= L[size - 1, size - 1];                        
            SearchResult();
            for (int i = 0; i < size; i++)
                Console.WriteLine("x{0} = {1}",i+1,res[i]);
            Check(matrixConst,fmembConst);
        }

        public static void Cholesky() {
            for (int i = 0; i < size; i++)
                fmemb[i] = 0;
            MatrixOnZero(matrix);
            MultipleMatrix();
            WriteMatrix(matrix, fmemb);
            MatrixOnZero(L1);            
            MatrixL();
            for(int i = 0; i<size;i++)
            Console.WriteLine("x{0} = {1}",i+1,res[i]);
            Check(matrix,fmemb);
        }
        public static void ReFactor(int step) {            
            
                for (int i = size-2; i >= step; i--)
                {
                    if (L[i, step] < L[i + 1, step])
                    {
                        double temp;
                        for (int j = 0; j < size; j++)
                        {
                            temp = L[i, j];
                            L[i, j] = L[i + 1, j];
                            L[i + 1, j] = temp;
                        }
                        temp = fmemb[i];
                        fmemb[i] = fmemb[i + 1];
                        fmemb[i + 1] = temp;
                    }

                }
            WriteMatrix(L, fmemb);
           
        }

        public static void Substraction(int step) {
            for (int i = step+1; i < size; i++)
            {
                for (int j = step; j < size; j++)
                {
                    L[i, j] -= L[step, j];
                }
                fmemb[i] -= fmemb[step];
            }
            WriteMatrix(L, fmemb);
        }

        public static void Division (int step) {
            for (int i = step; i < size; i++)
            {
                
                fmemb[i] /= L[i, step];
                for (int j = size-1; j >= 0; j--)
                {
                    L[i, j]/= L[i, step];
                }

                
            }
            WriteMatrix(L,fmemb);
        }
        public static void SearchResult() {

            for (int i = size - 1; i >= 0; i--) {
                double resuL1 = fmemb[i];
                for (int j = size - 1; j >= 0; j--) {
                    if (L[i, j] != 1)  
                        resuL1 -= L[i, j] * res[j];
                    else if (L[i, j] == 1) 
                        break;                    
                }
                res[i] = resuL1;
            }
        
        }
        
        public static void WriteMatrix(double[,] matx, double[] vec) {
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    Console.Write(matx[i, j] + "\t");
                }
                Console.WriteLine("||"+vec[i]);
                              
            }
            Console.WriteLine();
        }

        public static void MultipleMatrix() {
            
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    L[j, i] = matrixConst[i, j];
                    
                }
            }
          
            for (int i = 0; i < size; i++) {                
                for (int j = 0; j < size; j++) {                    
                    for (int p = 0; p < size; p++) {
                        matrix[i,j] += matrixConst[i, p] * L[p, j];
                    }
                }
            }
            
            for (int i = 0; i < size; i++)
            {
                fmemb[i] = 0;
                for (int j = 0; j < size; j++)
                {
                    fmemb[i] += L[i, j] * fmembConst[j];
                }
                                
            }
            

        }
        public static void MatrixL()
        {
            L[0, 0] = Math.Sqrt(matrix[0, 0]);

            for (int i = 1; i < size; i++)
            {
                double SummaJElementov = 0;
                for (int j = 1; j < size + 1; j++)
                {
                    if (j - 1 < i)
                    {
                        for (int k = 0; k < j - 1; k++)
                        {
                            SummaJElementov = SummaJElementov + L[i, k] * L[j - 1, k];
                        }
                        L[i, j - 1] = (matrix[i, j - 1] - SummaJElementov) / L[j - 1, j - 1];
                    }
                }
                double SumI_El = 0;
                for (int k = 0; k < i; k++)
                {
                    SumI_El +=L[i, k] * L[i, k];
                }
                L[i, i] = Math.Sqrt(matrix[i, i] - SumI_El);
            }

            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    L1[i, j] = L[j, i];
                }
            }

            
            for (int i = 0; i < size; i++)
            {
                double f = 0;
                for (int j = 0; j < i; j++)
                {
                    f +=L[i, j] * y[j];
                }
                y[i] = (fmemb[i] - f) / L[i, i];
                
            }

            
            for (int i = size - 1; i >= 0; i--)
            {
               double r = 0;
                for (int j = size - 1; j > i; j--)
                {
                    r +=L1[i, j] * res[j];
                }
                res[i] = (y[i] - r) / L1[i, i];                
            }

            

        }

        public static void Check(double[,] arr,double[] fm) {
            Console.WriteLine("\nПроверка\n");
            for (int i = 0; i < size; i++) {
                double result = 0;
                for (int j = 0; j < size; j++) {
                    result += res[j] * arr[i, j];
                }
                Console.WriteLine(result+" = "+fm[i]);
            }
        }
        public static void MatrixOnZero(double[,] arr) {
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++) {
                    arr[i, j] = 0;
                    
                    
                }
            }
            
        }
    }
}
