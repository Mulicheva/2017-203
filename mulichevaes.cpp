#include "mulichevaes.h"

/**
 * ����� ������
 */
void mulichevaes::lab1()
{
    double t;
    //������ ���
    for (int k = 0; k < N; k++)
    {
        t = A[k][k];

        //����� ��� �������� k'�� ������ �� a[k][k] �������, ������ ��� �� ������������ � �� ����� �� �� �������� �������� 1 ���
       // �������� ���������� ����������
        for (int j = 0; j < N; j++)
            A[k][j] = A[k][j] / t;
            b[k] =b[k]/t;

        for (int i = k + 1; i < N; i++)
        {
            t = A[i][k];
            //�������� �� ���� ����� ������� ���� k'�� �'�� ������ ����������� �� k'�� ������� ������
           // �� ������� ��������, ��� ��� ��� ���� � ���� �������� ����� ��������� � ������� ���������� ����������� ����������� ���
        
            for (int j = 0; j < N; j++)
            {
                A[i][j] =A[i][j]- A[k][j] * t;
            }
            b[i] =b[i]- b[k] * t;
        }
    }

    //������� ������������ ���� ������, ������ ����� ���������������, ����� ����� ��������� ������� �������� ��������� ������� x
    //������������ �������� ���
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
 * ����� ������ � ������� �������� ��������
 */
void mulichevaes::lab2()
{
//��� ��������� ����������� ������ ������. 
//���� ������� � ���, ��� �� ��� ������ ����, ������� ������� ���������, ��� ������� �� ����� �������� ���� � ����������� �������, 
//������ ���� ������������ �� ������ ����� �� ��������� ����������� �����. ������ ��� ��� ������� ��� ������� ��������� ������ ����������,
// �� � ������ ��� ����� ��������� ������� �� ����
    double t;
    
    for (int k = 0; k < N; k++)
    {
    //���� ����� ������� ������� ������� �� �'���� ����� �� ���� �������, �� ��������� �� k'��
        int maxel=k;
        for(int i=k+1;i<N;i++)
            if(abs(A[i][k]) > abs(A[maxel][k]))
                maxel=i;
     //�����, ������ ���� �������� k'�� ������� � ������� � ������������ ��������� �������
        for(int i=0;i<N;i++)
            std::swap(A[k][i],A[maxel][i]);
        std::swap(b[k],b[maxel]);


        //� ������ �� ������ �� ������ ������ �� ���� �� ������������ ����� 1

    //������ ���
   
        t = A[k][k];

 //����� ��� �������� k'�� ������ �� a[k][k] �������, ������ ��� �� ������������ � �� ����� �� �� �������� �������� 1 
//��� �������� ���������� ����������
        for (int j = 0; j < N; j++)
            A[k][j] = A[k][j] / t;
            b[k] =b[k]/t;

        for (int i = k + 1; i < N; i++)
        {
            t = A[i][k];
   //�������� �� ���� ����� ������� ���� k'�� �'�� ������ ����������� �� k'�� ������� ������,
// �� ������� ��������, ��� ��� ��� ���� � ���� �������� ����� ��������� � ������� ���������� ����������� ����������� ���
        
            for (int j = 0; j < N; j++)
            {
                A[i][j] =A[i][j]- A[k][j] * t;
            }
            b[i] =b[i]- b[k] * t;
        }
    }

    //������� ������������ ���� ������, ������ ����� ���������������, ����� ����� ��������� ������� �������� ��������� ������� x
    //������������ �������� ���
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
 * ����� ���������
 */
void mulichevaes::lab3()
{
	//����� ������� � ������������ � �������������? , ����� ��� ����������� � ���� A=L*LT, ��� LT ��� ����������������� ������� L
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
    //��������� D � L
    //L ����� ����������� �������, ��������� �������� �������� ������
    
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
    //������������� L
    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            LT[i][j] = L[j][i];
    //�������� LT �� D
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            LT[i][j] *= D[j];
         // ������ ������ ����� ������� Ax=B ����� �������� ������� �� ���� LT*y=B,L*x=y ����� �������  ������� ������ y, � ����� ��� ������� x
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
 * ����� ��������
 */


void mulichevaes::lab4()
{

//����� �������� - ������� ������ ������ ������.
//������� ���������������
    double *d, *c, *a; // �������, �������, ������ ���������
    d = new double[N];
    c = new double[N];
    a = new double[N];
    double m;

    a[0] = 0; 
    d[N-1] = 0;
    for (int i = 0; i < N; ++i)// ���������� ���������
    {
        if (i - 1 >= 0 && i - 1 < N) 
            d[i] = A[i-1][i]; //�������
        c[i] = A[i][i];      //�������
        if (i + 1 >= 0 && i + 1 < N) 
            a[i] = A[i+1][i]; //������
    }  
//���� ������ ��������
    for (int i = 1; i < N; i++) //����������� ����������� ������������
    {
        m = a[i]/c[i-1];
        c[i] = c[i] - m*d[i-1];
        b[i] = b[i] - m*b[i-1];
    }
//���� �������� �������� 
    x[N-1] = b[N-1]/c[N-1]; //����������� �������

    for (int i = N - 2; i >= 0; i--) //������� ��� ���������� ������
        x[i]=(b[i]-d[i]*x[i+1])/c[i];  

    delete[] b, c, d;
}



/**
 * ����� �����
 */
void mulichevaes::lab5()
{
    double *xold = new double[N];
    for (int i=0; i<N; i++)
    {
        x[i]=0; // �������������� ����� �������
    }
    double p=0.0;
    double eps=1e-20;
    int k=0;
    do
    {
        k++;
        p=0.0;
        for(int i=0; i<N; i++)
            xold[i]=x[i]; // ����� ������������ ���������� �������
        for(int i=0; i<N; i++)
        {
            double s=0; //��������� s, �� �� �� ���� ������������ ��������
            for(int j=0; j<i; j++)
                s += A[i][j] * xold[j];
            for(int j=i+1; j<N; j++)
                s += A[i][j] * xold[j];
            x[i]=(b[i] - s)/A[i][i]; // ����������� ����� ������� 
        }
        p= std::abs(xold[0]-x[0]);
        for(int i=0; i<N; i++)
        {
            if(std::abs(xold[i]-x[i]) > p) 
                p = std::abs(xold[i]-x[i]);//������������ ������� ����� ���������� �������� � �������.
        }
    } while(p >= eps);
    std::cout << "��c�� �������� : " << k << std::endl;

    delete [] xold;
}





/**
 * ����� �������
 */
void mulichevaes::lab6()
{
    double *xold = new double[N];
    for (int i=0; i<N; i++)
    {
        x[i]=0; // ��������� �����������
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
    std::cout << "��c�� �������� : " << k << std::endl;

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
  return std::string("�������� ���������");
}
