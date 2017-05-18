#include "mulichevaes.h"

/**
 * Метод Гаусса
 */
void mulichevaes::lab1()
{
    double t;
    //Прямой ход
    for (int k = 0; k < N; k++)
    {
        t = A[k][k];

        //Делим все элементы k'ой строки на a[k][k] элемент, потому что он диагональный и мы хотим на нём получить значение 1 для
       // удобства дальнейших вычислений
        for (int j = 0; j < N; j++)
            A[k][j] = A[k][j] / t;
            b[k] =b[k]/t;

        for (int i = k + 1; i < N; i++)
        {
            t = A[i][k];
            //Вычитаем из всех строк лежащих ниже k'ой к'ую строку помноженную на k'ый элемент строки
           // из которой вычитаем, что даёт нам ноль в этом элементе после вычитания и матрица постепенно приобретает треугольный вид
        
            for (int j = 0; j < N; j++)
            {
                A[i][j] =A[i][j]- A[k][j] * t;
            }
            b[i] =b[i]- b[k] * t;
        }
    }

    //Матрица треугольного вида готова, теперь можем последовательно, снизу вверх вычислять искомые значения элементов матрицы x
    //Осуществляем обратный ход
    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            t = A[i][k];
            
            for (int j = 0; j < N; j++)
                A[i][j] =A[i][j]- A[k][j] * t;
            b[i] =b[i] - b[k] * t;
        }
    }

    for(int i=0; i<N; i++)
        x[i]=b[i];
    
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void mulichevaes::lab2()
{
//Это некоторая оптимизация метода Гаусса. 
//Суть состоит в том, что бы при каждом шаге, элемент главной диагонали, под которым мы будем получать нули в нижележащих строках, 
//должен быть максимальным из других таких же элементов нижележащих строк. Потому как это поможет при делении уменьшить ошибку округления,
// ну и вообще там может оказаться деление на ноль
    double t;
    
    for (int k = 0; k < N; k++)
    {
    //Ищем самый большой элемент стоящий на к'атом месте по всем строкам, по умолчанию он k'ый
        int maxel=k;
        for(int i=k+1;i<N;i++)
            if(abs(A[i][k]) > abs(A[maxel][k]))
                maxel=i;
     //Нашли, теперь надо поменять k'ую строчку и строчку с максимальным элементом местами
        for(int i=0;i<N;i++)
            std::swap(A[k][i],A[maxel][i]);
        std::swap(b[k],b[maxel]);


        //А дальше всё делаем по методу Гаусса по коду из лабораторной номер 1

    //Прямой ход
   
        t = A[k][k];

 //Делим все элементы k'ой строки на a[k][k] элемент, потому что он диагональный и мы хотим на нём получить значение 1 
//для удобства дальнейших вычислений
        for (int j = 0; j < N; j++)
            A[k][j] = A[k][j] / t;
            b[k] =b[k]/t;

        for (int i = k + 1; i < N; i++)
        {
            t = A[i][k];
   //Вычитаем из всех строк лежащих ниже k'ой к'ую строку помноженную на k'ый элемент строки,
// из которой вычитаем, что даёт нам ноль в этом элементе после вычитания и матрица постепенно приобретает треугольный вид
        
            for (int j = 0; j < N; j++)
            {
                A[i][j] =A[i][j]- A[k][j] * t;
            }
            b[i] =b[i]- b[k] * t;
        }
    }

    //Матрица треугольного вида готова, теперь можем последовательно, снизу вверх вычислять искомые значения элементов матрицы x
    //Осуществляем обратный ход
    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            t = A[i][k];
            
            for (int j = 0; j < N; j++)
                A[i][j] =A[i][j]- A[k][j] * t;
            b[i] =b[i] - b[k] * t;
        }
    }

    for(int i=0; i<N; i++)
        x[i]=b[i];
}


/**
 * Метод Холецкого
 */
void mulichevaes::lab3()
{
	//Пусть матрица А симметричная и положительная? , тогда она представима в виде A=L*LT, где LT это транспонированная матрица L
    double **L = new double*[N];
    double **LT = new double*[N];
    double *D = new double[N];
    for(int i=0; i<N; i++)
    {
       L[i] = new double[N];
       LT[i] = new double[N];
       D[i]=0;
       for(int j=0; j<N; j++)
       {
           L[i][j]=0;
       }
    }
    //вычисляем D и L
    //L нижне треугольная матрица, находится согласно формулам метода
    
    for(int i=0; i<N; i++)
    {
        double isum = 0;
        for(int k=0; k<i; k++)
                isum += std::abs(L[k][i]) * std::abs(L[k][i]) * D[k];

        if((A[i][i] - isum) > 0)
            D[i] = 1;
        else
            D[i] = -1;

        L[i][i] = sqrt(std::abs(A[i][i] - isum));

        for(int j=i+1; j<N; j++)
        {
            double sum=0;
            for(int k=0; k<i; k++)
                sum += L[k][i] * L[k][j] * D[k];
            L[i][j] = (A[i][j] - sum) / (L[i][i] * D[i]);
        }
    }
    //транспонируем L
    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            LT[i][j] = L[j][i];
    //умножаем LT на D
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            LT[i][j] *= D[j];
         // Теперь вместо одной матрицы Ax=B можем получить систему из двух LT*y=B,L*x=y найдём решение  выразив сперва y, а потом уже искомый x
    //LT*D * y = b
    double *y = new double[N];
    for(int i=0; i<N; i++)
    {
        double s=0;
        for(int j=0; j<i; j++)
            s += y[j] * LT[i][j];
        y[i]=(b[i]-s)/LT[i][i];
    }
    //L * x = y
    for(int i=N-1; i>=0; i--)
    {
        double s=0;
        for(int j=i+1; j<N; j++)
            s += x[j] * L[i][j];
        x[i]=(y[i] - s)/L[i][i];
    }

    for(int i=0; i<N; i++)
    {
        delete [] L[i];
        delete [] LT[i];
    }
    delete [] L;
    delete [] LT;
    delete [] D;
    delete [] y;
}




/**
 * Метод Прогонки
 */


void mulichevaes::lab4()
{

//Метод прогонки - частный случай метода Гаусса.
//Матрица трёхдиагональная
    double *d, *c, *a; // верхняя, главная, нижняя диагонали
    d = new double[N];
    c = new double[N];
    a = new double[N];
    double m;

    a[0] = 0; 
    d[N-1] = 0;
    for (int i = 0; i < N; ++i)// выписываем диагонали
    {
        if (i - 1 >= 0 && i - 1 < N) 
            d[i] = A[i-1][i]; //верхняя
        c[i] = A[i][i];      //главная
        if (i + 1 >= 0 && i + 1 < N) 
            a[i] = A[i+1][i]; //нижняя
    }  
//Этап прямой прогонки
    for (int i = 1; i < N; i++) //Вычисляются прогоночные коэффициенты
    {
        m = a[i]/c[i-1];
        c[i] = c[i] - m*d[i-1];
        b[i] = b[i] - m*b[i-1];
    }
//Этап обратной прогонки 
    x[N-1] = b[N-1]/c[N-1]; //вычисляется решение

    for (int i = N - 2; i >= 0; i--) //решение для дальнейших частей
        x[i]=(b[i]-d[i]*x[i+1])/c[i];  

    delete[] b, c, d;
}



/**
 * Метод Якоби
 */
void mulichevaes::lab5()
{
    double *xold = new double[N];
    for (int i=0; i<N; i++)
    {
        x[i]=0; // первоначальное новое решение
    }
    double p=0.0;
    double eps=1e-20;
    int k=0;
    do
    {
        k++;
        p=0.0;
        for(int i=0; i<N; i++)
            xold[i]=x[i]; // здесь записывается предыдущее решение
        for(int i=0; i<N; i++)
        {
            double s=0; //вычисляем s, но мы не берём диагональные элементы
            for(int j=0; j<i; j++)
                s += A[i][j] * xold[j];
            for(int j=i+1; j<N; j++)
                s += A[i][j] * xold[j];
            x[i]=(b[i] - s)/A[i][i]; // вычисляется новое решение 
        }
        p= std::abs(xold[0]-x[0]);
        for(int i=0; i<N; i++)
        {
            if(std::abs(xold[i]-x[i]) > p) 
                p = std::abs(xold[i]-x[i]);//максимальная разница между предыдущим решением и текущим.
        }
    } while(p >= eps);
    std::cout << "Чиcло итераций : " << k << std::endl;

    delete [] xold;
}





/**
 * Метод Зейделя
 */
void mulichevaes::lab6()
{
    double *xold = new double[N];
    for (int i=0; i<N; i++)
    {
        x[i]=0; // начальное приближение
    }
    double p=0.0;
    double eps=1e-20;
    int k=0;
    do
    {
        k++;
        p=0.0;
        for(int i=0; i<N; i++)
            xold[i]=x[i];
        for(int i=0; i<N; i++)
        {
            double s=0;
            for(int j=0; j<i; j++)
                s += A[i][j] * x[j];
            for(int j=i+1; j<N; j++)
                s += A[i][j] * xold[j];
            x[i]=(b[i] - s)/A[i][i];
        }
        for(int i=0; i<N; i++)
        {
            if(std::abs(xold[i]-x[i]) > p)
                p = std::abs(xold[i]-x[i]);
        }

    } while(p >= eps);
    std::cout << "Чиcло итераций : " << k << std::endl;

    delete [] xold;
}




/**
 * ???? ?? ??????????? ???????
 */
void mulichevaes::lab7()
{

}
void mulichevaes::lab8()
{

}

std::string mulichevaes::get_name()
{        
  return std::string("Муличева Екатерина");
}
