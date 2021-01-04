#include <iostream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;

// The number of Long numbers rhat was created
int NUM_COUNT = 0;

// Struct
struct Long
{
    string hex_str;
    long long big_num[199] = { 0 };
    int len = 0;
    int bin[9999] = { 0 };
};


// using in Input() function
int to_int(char n)
{
    switch (n)
    {
    case '0':
        return 0;
        break;
    case '1':
        return 1;
        break;
    case '2':
        return 2;
        break;
    case '3':
        return 3;
        break;
    case '4':
        return 4;
        break;
    case '5':
        return 5;
        break;
    case '6':
        return 6;
        break;
    case '7':
        return 7;
        break;
    case '8':
        return 8;
        break;
    case '9':
        return 9;
        break;
    case 'A':
        return 10;
        break;
    case 'B':
        return 11;
        break;
    case 'C':
        return 12;
        break;
    case 'D':
        return 13;
        break;        
    case 'E':
        return 14;
        break;        
    case 'F':
        return 15;
        break;        
    
    default:
        break;
    }
}


// Intup function
void Input(Long &N)
{
    if(N.hex_str.length() == 0)
    {
        cout << "Input number A" << ++NUM_COUNT << " in hex:" << endl;
        cin >> N.hex_str;
    }
    N.len = 0;
    int size = N.hex_str.length();
    int hex_coef = 0;
    int j = 0; //counter for elem of big_num
    int deg = 0;
    for (int i = size - 1; i >= 0; i--)
    {
        int symbol = N.hex_str[i];
        if((symbol >= 48 && symbol <= 57) || (symbol >= 65 && symbol <= 70))
        {
            if (j % 8 == 0 && j != 0)
            {
                N.len++;
                N.big_num[N.len] = 0;
            }
            hex_coef = to_int(N.hex_str[i]);
            deg = (size - i - 1) % 8;
            N.big_num[N.len] = N.big_num[N.len] + hex_coef * pow(16,deg);   
            j++;
        }
        else
        {
            cout << "Wrong input!!!" << endl;
            exit(0);
        }   
    }
    N.len++;
}


// Output function
void Output(Long N)
{
    N.hex_str = "";
    long long tmp = 0;
    char ch = ' ';
    for (int i = 0; i < N.len; i++)
    {
        for (int j = 7; j >= 0; j--)
        {
            tmp = N.big_num[i] % 16;
            if( tmp >= 0 && tmp <= 9)
            {
                ch = tmp + 48;
                N.hex_str = ch + N.hex_str;
            }
            if (tmp >= 10 && tmp <= 15)
            {
                ch = tmp + 55;
                N.hex_str = ch + N.hex_str;
            }
            N.big_num[i] = N.big_num[i] / 16;
        }
    }
    cout << "Result = ";
    int k = 0;
    if(N.hex_str[0] == '0' && N.hex_str.length() > 1)
    {
        while (N.hex_str[k] == '0')
        {
            k++;
        }
        int tmp_size = N.hex_str.length();
        int * hex_tmp = new int[tmp_size];

        for (int i = k; i < tmp_size; i++)
        {
            hex_tmp[i - k] = N.hex_str[i];
        }
        N.hex_str = "";
        for (int i = 0; i < tmp_size - k; i++)
        {
            N.hex_str.push_back(hex_tmp[i]);
        }
    }        
        
    for (int i = 0; i < N.hex_str.length(); i++)
    {
        cout << N.hex_str[i];
    }
    cout << endl;
}




string Output_str(Long N)
{
    N.hex_str = "";
    long long tmp = 0;
    char ch = ' ';
    for (int i = 0; i < N.len; i++)
    {
        for (int j = 7; j >= 0; j--)
        {
            tmp = N.big_num[i] % 16;
            if( tmp >= 0 && tmp <= 9)
            {
                ch = tmp + 48;
                N.hex_str = ch + N.hex_str;
            }
            if (tmp >= 10 && tmp <= 15)
            {
                ch = tmp + 55;
                N.hex_str = ch + N.hex_str;
            }
            N.big_num[i] = N.big_num[i] / 16;
        }
    }
    cout << "Result = ";
    int k = 0;
    if(N.hex_str[0] == '0' && N.hex_str.length() > 1)
    {
        while (N.hex_str[k] == '0')
        {
            k++;
        }
        int tmp_size = N.hex_str.length();
        int * hex_tmp = new int[tmp_size];

        for (int i = k; i < tmp_size; i++)
        {
            hex_tmp[i - k] = N.hex_str[i];
        }
        N.hex_str = "";
        for (int i = 0; i < tmp_size - k; i++)
        {
            N.hex_str.push_back(hex_tmp[i]);
        }
    }

    return N.hex_str;            
    cout << endl;
}

// Add function
Long Add(Long A, Long B)
{
    Long Res;
    unsigned long long res_coef = 0;
    int carry = 0;
    long long base = 4294967296;
    unsigned long long tmp = 0;
    Res.len = max(A.len,B.len);
    for (int i = 0; i < Res.len; i++)
    {

        res_coef = A.big_num[i] + B.big_num[i] + carry;
        Res.big_num[i] = res_coef % base;
        tmp = res_coef / base;
        carry = tmp;

        /*
        res_coef = A.big_num[i] + B.big_num[i] + carry;
        Res.big_num[i] = res_coef & (b - 1);
        tmp = (res_coef >> 32);
        carry = tmp;*/
    }
    
    if(carry)
    {
        Res.big_num[Res.len] = carry;
        Res.len++;
    }

    return Res;
}

// Comparison A Long and B Long
int Compare(Long A, Long B)
{
    if(A.len > B.len) 
        return 1;
    if(A.len < B.len)
        return -1;
    for (int i = A.len - 1; i >= 0; i--)
    {
        if(A.big_num[i] > B.big_num[i])
        {
            return 1;
        }
        if (A.big_num[i] < B.big_num[i])
        {
            return -1;
        }
    }
    return 0;
}

//Subtraction function
Long Sub(Long A, Long B)
{
    Long Res;
    long long base = 1;
    int cmp = Compare(A,B);
    long long borrow = 0;
    long long tmp = 0;
    base = (base << 32);
    if(cmp == 1)
    {
        for (int i = 0; i < A.len; i++)
        {
            tmp = A.big_num[i] - B.big_num[i] - borrow;
            if(tmp >= 0)
            {
                Res.big_num[i] = tmp;
                Res.len++;
                borrow = 0;
            }
            else
            {
                Res.big_num[i] = base + tmp;
                Res.len++;  
                borrow = 1;
            }
            
        }
        while (Res.big_num[Res.len - 1] == 0)
        {
            Res.len--;
        }
        
        return Res;
    }
    if(cmp == -1)
    {
        cout << "Negative number" << endl;
        return Res;
    }
    if(cmp == 0)
    {
        Res.big_num[0] = 0;
        Res.len = 1;
        Res.hex_str = "0";
        return Res;
    }
    return Res;
}

// Multipliaction Of Long A on Long B[i] function
// Using in Mult(Long,Long) function
Long MultOneDig(Long A, long long b)
{
    Long B;
    long long base = 1;
    unsigned long long carry = 0;
    unsigned long long tmp = 0;
    base = (base << 32);
    B.len = 0;
    *B.big_num = { 0 };
    for (int i = 0; i < A.len; i++)
     {
        tmp = A.big_num[i] * b + carry;
		B.big_num[i] = tmp % base;
		B.len++;
        carry = tmp / base;
		
     }

    //     tmp = A.big_num[i] * b + carry;
    //     B.big_num[i] = tmp & (base - 1);
    //     B.len++;
    //     carry = (tmp >> 32);
    // 

    if(carry)
    {
        B.big_num[B.len] = carry;
        B.len++;
    }
    
    return B;
}

//Shifting Long Num elements
void ShiftDigToHigh(Long &A, int pos)
{
    if(pos > 0)
    {
        Long C;
        for (int i = 0; i < A.len; i++)
        {
            C.big_num[i] = A.big_num[i];
            C.len++;
        }
        for (int i = 0; i < C.len; i++)
        {
            A.big_num[i + pos] = C.big_num[i];
        }

        for (int i = 0; i < pos; i++)
        {
            A.big_num[i] = 0;
        }
        A.len += pos;      
    }
}

// Multiplication func
Long Mult(Long A, Long B)
{
    Long TmpMult;
    Long Res;
    
    int length = max(A.len,B.len);
    for (int i = 0; i < length; i++)
    {
        TmpMult = MultOneDig(A, B.big_num[i]);
        ShiftDigToHigh(TmpMult, i);
        Res = Add(Res, TmpMult);
    }

    return Res;
}


// Len of Bits
int BitLen(Long N)
{
    
    // int tmp = 0;
    int len = 0;
    len = log2(N.big_num[N.len - 1]) + 1;
    len += (N.len - 1) * 32;
    return len;
}

Long ShiftBitsToHigh(Long N, int pos)
{
   // Long Tmp;
    for (int i = 0; i < pos; i++)
    {
        N = MultOneDig(N, 2);
    } 
    return N;
}

Long SetBits(Long N, int pos)
{
    while (N.bin[pos] == 1)
    {
        N.bin[pos] = 0;
        pos++;
    }

    N.bin[pos] = 1;
    
    if(N.len <= pos)
    {
        N.len = pos + 1;
    }
    
    return N;
}

// binary to hex
void bin_to_hex(Long &N)
{
    //char hex_ch = '\0';
    unsigned long long *tmp = new unsigned long long[N.len];
    int i = 0;
    int k = 0;
    unsigned long long elem = 0;
    unsigned long long tmp_el = 0;
    unsigned long long sum = 0;

    while (i < N.len)
    {
        for (int j = i; j < i + 4; j++)
        {
            elem = N.bin[j];
            tmp_el = pow(2, j % 4) * elem;
            sum += tmp_el;
            tmp[k] = sum;
        }
        sum = 0;
        k++;
        i += 4;
    }

    //int j = 0;
    unsigned long long hex_tmp = 0;
    for (int i = k - 1; i >= 0; i--)
    {
       // j++;
        hex_tmp = tmp[i];
        switch (hex_tmp)
        {
       case 0:
//            N.hex_str[j] = '0';
            N.hex_str = N.hex_str + '0';
            break;
       case 1:
         //   N.hex_str[j] = '1';
            N.hex_str = N.hex_str + '1';
            break;
        case 2:
           // N.hex_str[j] = '2';
            N.hex_str = N.hex_str + '2';
            break;
        case 3:
           // N.hex_str[j] = '3';
            N.hex_str = N.hex_str + '3';
            break;
        case 4:
          //  N.hex_str[j] = '4';
            N.hex_str = N.hex_str + '4';
            break;
        case 5:
          //  N.hex_str[j] = '5';
            N.hex_str = N.hex_str + '5';
            break;
        case 6:
         //   N.hex_str[j] = '6';
            N.hex_str = N.hex_str + '6';
            break;
        case 7:
        //    N.hex_str[j] = '7';
            N.hex_str = N.hex_str + '7';
            break;
        case 8:
        //    N.hex_str[j] = '8';
            N.hex_str = N.hex_str + '8';
            break;
        case 9:
        //    N.hex_str[j] = '9';
            N.hex_str = N.hex_str + '9';
            break;
        case 10:
       //     N.hex_str[j] = 'A';
            N.hex_str = N.hex_str + 'A';
            break;
        case 11:
      //      N.hex_str[j] = 'B';
            N.hex_str = N.hex_str + 'B';
            break;
        case 12:
     //       N.hex_str[j] = 'C';
            N.hex_str = N.hex_str + 'C';
            break;
        case 13:
      //      N.hex_str[j] = 'D';
            N.hex_str = N.hex_str + 'D';
            break;
        case 14:
      //      N.hex_str[j] = 'E';
            N.hex_str = N.hex_str + 'E';
            break;
        case 15:
     //       N.hex_str[j] = 'F';
            N.hex_str = N.hex_str + 'F';
            break;
        
        default:
            break;
        }
    }
    //N.len = j;
    cout << endl;    

}


// Div func
void Div(Long A,Long B,Long &R, Long &Q)
{
    if(B.big_num[B.len - 1] != 0)
    {
        Long Tmp;
        Long Tmp1;
        int B_bit_len = BitLen(B);
        int R_bit_len;
        R = A;
        int cmp_R_B = Compare(R, B);
        int cmp_R_Tmp;
        while (cmp_R_B >= 0)
        {
            R_bit_len = BitLen(R);
            Tmp = ShiftBitsToHigh(B, R_bit_len - B_bit_len);
            cmp_R_Tmp = Compare(R, Tmp);
            while (cmp_R_Tmp == -1)
            {
                R_bit_len--;
                Tmp = ShiftBitsToHigh(B, R_bit_len - B_bit_len);
                cmp_R_Tmp = Compare(R, Tmp);
            }
            R = Sub(R, Tmp);
            Q = SetBits(Q, R_bit_len - B_bit_len);

            cmp_R_B = Compare(R, B);
        }
        

    }

}

// write bin array

Long Write_bin(Long N)
{
    int c;
    N.len = 0;
    int coef = 0;
    int len_N_hex = N.hex_str.length();
    for (int i = len_N_hex - 1; i >= 0; i--)
    {   
        coef = N.hex_str[i];
        
        if (coef >= '0' && coef <= '9')
        {
            c = coef - 48;       
        }
        if (coef >= 'A' && coef <= 'F')
        {
            c = coef - 55;
        }
        
        for (int i = 0; i < 3; i++)
        {
            N.bin[N.len] = c % 2;
            N.len++;
            c /= 2;    
        }
        
        N.bin[N.len] = c;
        N.len++;
    }

    while (N.bin[N.len - 1] == 0 && N.len > 1)
    {
        N.len--;
    }
    return N;
    
}




// Long Power

Long Pow(Long A, Long B)
{
    Long C;
    if (B.bin[B.len - 1] == 0)
    {
        C.len = 1;
        C.big_num[0] = 1;
    }
    else
    {
        C = A;
        for (int i = B.len - 2; i >= 0; i--)
        {
            C = Mult(C,C);

            if (B.bin[i] == 1)
            {
                C = Mult(C,A);
            }
            
        }
        
    }

    while (C.big_num[C.len - 1] == 0)
    {
        C.len--;
    }
    return C; 
    
}



void Distribution(Long A, Long B, Long C)
{
    Long N1;
    Long N2;
    Long N3;
    Long N4;
    Long R1;
    Long R2;
    Long R3;
    N1 = Add(A,B);
    R1 = Mult(N1,C);
    cout << "(A + B) * C :" <<endl; 
    Output(R1);

    R2 = Mult(C,N1);
    cout << " C * (A + B) :" <<endl; 
    Output(R2);

    N3 = Mult(A,C);
    N4 = Mult(B,C);
    R3 = Add(N3,N4);
    cout << "A * C + B * C :" << endl;
    Output(R3);
}


void Constant(Long A)
{
    int n = 150;
    Long B;
    B = Add(A,A);
    cout<< "\n iteration: 1 :" <<endl;
    Output(B);
    
    for (int i = 0; i < n; i++)
    {
        cout<< "\n iteration: " << i + 2 <<" :"<<endl;
        B = Add(A,B);
        Output(B);
        
    }

    cout<<endl;
    Long Ten;
    Long D;
    Input(Ten);
    D = Mult(A,Ten);
    Output(A);
    Output(D);
}






// LAB2:

bool Parity(Long N)
{
    int p = N.big_num[0] & 1;
    if(p == 0)
        return true;
    return false;
}

Long DivOnTwo(Long N)
{
    unsigned long int base = 4294967296;
    for (int i = 0; i < N.len; i++)
    {
        N.big_num[i] = (N.big_num[i] >> 1) + (N.big_num[i+1] << 31) & (base - 1);
    }

    while (N.big_num[N.len - 1] == 0)
    {
        N.len--;
    }
    return N;
    
}


Long GCD(Long A, Long B)
{
    Long D;
    Long Tmp;
    D.big_num[0] = 1;
    D.len = 1;

    int cmp = 0;
    
    while (Parity(A) && Parity(B))
    {
        A = DivOnTwo(A);
        B = DivOnTwo(B);
        D = MultOneDig(D, 2);
    }
    while (B.big_num[B.len - 1] != 0)
    {
        while (Parity(B))
        {
            B = DivOnTwo(B);
        }
        cmp = Compare(A,B);
        if (cmp >= 0)
        {
            Tmp = A;
            A = B;
            B = Sub(Tmp,B);
        }
        else
        {
            B = Sub(B,A);
        }
    }

    D = Mult(A,D);
    while (D.big_num[D.len - 1] == 0)
    {
        D.len--;
    }
    return D;
}


// LCM(a,b) using GCD:
// Input:
//     | a,b
// Output:
//     | (a * b) / GCD(a,b)

Long LCM(Long A, Long B)
{
    Long Lcm;
    Long Gcd = GCD(A,B);
    if(Gcd.big_num[0] == 1)
    {
        return Mult(A,B);
    }
    else
    {
        Long Q,R,R_tmp;
        Div(A,Gcd,R,Q);
        bin_to_hex(Q);
        Input(Q);
        Div(B,Gcd,R_tmp, R);
        bin_to_hex(R);
        Input(R);
        R_tmp = Mult(R,Q);
        Lcm = Mult(R_tmp,Gcd);
    }
    return Lcm;
}


// Long LCM(Long A, Long B)
// {
//     Long 

// }


Long AddMod(Long A, Long B, Long N)
{
    Long Res,Q, Sum;
    Sum = Add(A,B);
    Div(Sum, N, Res, Q);
    return Res;
}

Long SubMod(Long A, Long B, Long N)
{
    Long Res, Q, Subt, C;
    int cmp_AB = Compare(A,B);
    if(cmp_AB == -1 )
    {
        C = A;
        A = B;
        B = C; 
    }
    Subt = Sub(A,B);
    if(Compare(Subt,N) >= 0)
    {
        Div(Subt, N, Res, Q);
        return Res;
    }
    return Subt;

}



Long KillLastDig(Long A, int shift)
{
    Long Tmp;
    Tmp.len = 0;
    for (int i = 0; i < A.len-shift; i++)
    {
        Tmp.big_num[i] = A.big_num[i+shift];
        Tmp.len++;
    }
    
    while (Tmp.big_num[Tmp.len-1]==0 && Tmp.len > 1)
    {
        Tmp.len--;
    }
    
    return Tmp;
}


Long MultModBase(Long A, int BasePow)
{
    Long Tmp;
    for (int i = 0; i < BasePow; i++)
    {
        Tmp.big_num[i] = A.big_num[i];
        Tmp.len++;
    }

    while (Tmp.big_num[Tmp.len-1] == 0 && Tmp.len > 1)
    {
        Tmp.len--;
    }
    return Tmp;
    
    
}

Long BarretReduc(Long X, Long N, Long M)
{
    Long Q1,Q2,Q3, R,R1,R2, TmpMult;
    Q1 = KillLastDig(X, N.len-1);
    Q2 = Mult(Q1,M);
    Q3 = KillLastDig(Q2, N.len+1);
    
    R1 = MultModBase(X, N.len+1);
    TmpMult = Mult(Q3,N);
    R2 = MultModBase(TmpMult,N.len+1);

    
    R = Sub(R1,R2);
    while (Compare(R,N) >= 0)
    {
        R = Sub(R,N);
    }
    return R;
}

void MultModProc(Long A,Long B, Long N)
{
    Long Res, X, M, Base, Tmp;
    X = Mult(A,B);
    X.len = 2 * N.len;
    Base.len = 2;
    Base.big_num[1] = 1;
    Base.big_num[0] = 0;
    ShiftDigToHigh(Base, X.len - 1);
    Div(Base, N, Tmp, M);
    bin_to_hex(M);
    M.len = 0;
    Input(M);
    double start = clock();
    Res = BarretReduc(X, N, M);
    double end = clock();
    cout << "BArret time = " << (end - start) / 1000 << endl;
    Output(Res);
}


int main()
{
     unsigned int start_time = 0; // начальное время
    // здесь должен быть фрагмент кода, время выполнения которого нужно измерить
     unsigned int end_time = 0; // конечное время
     double search_time = 0; // искомое время    


    Long A;
    Long B;
    Long N;
    Long R;
    Long Add_res;
    Long Sub_res;
    Long Mult_res;
    Long Div_res;
    Long Pow_res;

    Long AddMod_res;
    Long SubMod_res; 
    Long GCD_res;
    Long LCM_res;

//Power
    
    
    A.hex_str="94EDE1A444B9738ADF06CDB40DCAFA87B25A8BECA2D2262A53D8431A119405F0CBEFB83D2AD547CCE3AE74A8EC74A313C8BED20D4349D9EFBA356FE6E8AD89E2";
    B.hex_str="5BCC0B222EE17877C9EB60FA91632BC7A6E29D80F02CD3FE16B5C2A2231B43DB2B2D12F21B293AAF49FE1165CB7A21D12D6ACEC225285544B36BABD3F8B4DD8D";
    N.hex_str="8E06E4DFFB37B57A66ECC52CF2D7D888C49C2794E6FB944C4183A128203932FEBEA4B6E62B2EBDAD4FF0B80DBEDC8439D31280D13E7E523596D92861F6A89E81";

    Input(A);
    Input(B);
    Input(N);


    // start_time = clock();
    // MultModProc(A,B,N);
    // end_time = clock();

    // // Output(AddMod_res);
    // search_time = end_time - start_time;
    // cout << ".time = " << search_time/1000 <<" sec"<< endl;

    // Write_bin(A);
    // Write_bin(B);
    // start_time = clock();
    // Div(A,B,R,Div_res);
    // end_time = clock();

    // bin_to_hex(Div_res);
    // cout << Div_res.hex_str << endl;
    // search_time = end_time - start_time;
    // cout << "Div.time = " << search_time/1000 <<" sec"<< endl;
    

    
    
    

    // start_time = clock();
    // LCM_res = LCM(A,B);
    // end_time = clock();
    // search_time = end_time - start_time;
    // cout << "LCM: " << endl;
    // Output(LCM_res);
    // cout << "LCM.time = " << search_time / 1000 <<" milisecond"<< endl;
    

    // start_time = clock();
    // GCD_res = GCD(A,B);
    // end_time = clock();
    // search_time = end_time - start_time;
    // cout << "GCD:" <<endl;
    // Output(GCD_res);
    // cout << "GCD.time = " << search_time / 1000 <<" milisecond"<< endl;
    
    // start_time = clock();
    // AddMod_res = AddMod(A,B,N);
    // end_time = clock();
    // search_time = end_time - start_time;
    // cout << "AddMod :" <<endl;
    // Output(AddMod_res);
    // cout << "AddMod.time = " << search_time / 1000<<" milisecond"<< endl;
    
    // start_time = clock();
    // SubMod_res = SubMod(A,B,N);
    // end_time = clock();
    // search_time = end_time - start_time;
    // cout << "SubMod:"<<endl;
    // Output(SubMod_res);
    // cout << "SubMod.time = " << search_time / 1000 <<" milisecond"<< endl;
    
    start_time = clock();
    MultModProc(A,B,N);
    end_time = clock();
    search_time = end_time - start_time;
    cout << "\nMultModProc.time = " << search_time <<" milisecond"<< endl;


//     start_time = clock();
//     Add_res = Add(A,B);
//     end_time = clock();
//     search_time = end_time - start_time;
//     cout << "Add.time = " << search_time <<" milisecond"<< endl;
    
//     start_time = clock();
//     Sub_res = Sub(A,B);
//     end_time = clock();
//     search_time = end_time - start_time;
//     cout << "Sub.time = " << search_time <<" milisecond"<< endl;
    
//     start_time = clock();
//     Mult_res = Mult(A,B);
//     end_time = clock();
//     search_time = end_time - start_time;
//     cout << "Mult.time = " << search_time <<" milisecond"<< endl;
// //Power   
//     Long X,Y;
//     X.hex_str = "ABC13";
//     Y.hex_str = "AB";
//     Input(X);
//     Input(Y);
//     X = Write_bin(X);
//     Y = Write_bin(Y);

//     start_time = clock();
//     Pow_res = Pow(A,B);
//     end_time = clock();
//     search_time = end_time - start_time;
//     cout << "Power.time = " << search_time <<" milisecond"<< endl;    
    
    
    
    
    


// //Add
//     cout << "Add(A,B): " << endl;
//     Add_res = Add(A,B);
//     Output(Add_res);
// //Sub
//     cout << "Sub(A,B): " << endl;
//     Sub_res = Sub(A,B);
//     Output(Sub_res);
// //Mult
//     cout << "Mult(A,B): " << endl;
//     Mult_res = Mult(A,B);
//     Output(Mult_res);
// //Div
//     Write_bin(A);
//     Write_bin(B);
//     Div(A,B,R,Div_res);
//     bin_to_hex(Div_res);
//     cout << "Div_res = " << Div_res.hex_str << endl; 



// //GCD
//     cout << "GCD(A,B) :" << endl;
//     Output(GCD_res);
// //LCM
//     cout << "LCM(A,B) :" << endl;
//     Output(LCM_res);
// //AddMod
//     cout << "A+B modN :" << endl;
//     Output(AddMod_res);
// //SubMod
//     cout << "A-B modN :" << endl;
//     Output(SubMod_res);
// //MultMod
//     cout << "A*B modN :" << endl;
//     MultModProc(A,B,N);

   
    return 0;
}

