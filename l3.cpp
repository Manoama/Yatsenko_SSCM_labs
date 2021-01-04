#include <iostream>
#include<string>
#include <cmath>
#include <ctime>

using namespace std;

struct Long
{
    string hex_str;
    int deg=0;
    long long coef[33] = { 0 };
    int max_bit=0;
    int bin[1200] = { 0 };
};

Long Gen;
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

void Input(Long &N)
{
    if(N.hex_str.length() == 0)
    {
        cout << "Input number A in hex:" << endl;
        cin >> N.hex_str;
    }
    N.deg = 0;
    int size = N.hex_str.length();
    int hex_coef = 0;
    int j = 0; //counter for elem of coef
    int deg = 0;
    for (int i = size - 1; i >= 0; i--)
    {
        int symbol = N.hex_str[i];
        if((symbol >= 48 && symbol <= 57) || (symbol >= 65 && symbol <= 70))
        {
            if (j % 8 == 0 && j != 0)
            {
                N.deg++;
                N.coef[N.deg] = 0;
            }
            hex_coef = to_int(N.hex_str[i]);
            deg = (size - i - 1) % 8;
            N.coef[N.deg] = N.coef[N.deg] + hex_coef * pow(16,deg);   
            j++;
        }
        else
        {
            cout << "Wrong input!!!" << endl;
            exit(0);
        }   
    }
    N.deg++;
}

// Output function
void Output(Long N)
{
    N.hex_str = "";
    long long tmp = 0;
    char ch = ' ';
    for (int i = 0; i < N.deg; i++)
    {
        for (int j = 7; j >= 0; j--)
        {
            tmp = N.coef[i] % 16;
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
            N.coef[i] = N.coef[i] / 16;
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

// Output function
void Output_str(Long &N)
{
    N.hex_str = "";
    long long tmp = 0;
    char ch = ' ';
    for (int i = 0; i < N.deg; i++)
    {
        for (int j = 7; j >= 0; j--)
        {
            tmp = N.coef[i] % 16;
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
            N.coef[i] = N.coef[i] / 16;
        }
    }
    //cout << "Result = ";
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
        
    // for (int i = 0; i < N.hex_str.length(); i++)
    // {
    //     cout << N.hex_str[i];
    // }
    // cout << endl;
}

void Write_bin(Long &N)
{
    int c;
    N.max_bit = 0;
    int coef = 0;
    int len_N_hex = N.hex_str.length();
    for (int i = len_N_hex - 1; i >= 0; i--)
    {   
        coef = N.hex_str[i];
        // make str_to_int funct
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
            N.bin[N.max_bit] = c % 2;
            N.max_bit++;
            c /= 2;    
        }
        
        N.bin[N.max_bit] = c;
        N.max_bit++;
    }

    while (N.bin[N.max_bit - 1] == 0 && N.max_bit > 1)
    {
        N.max_bit--;
    }   
}

void bin_to_hex(Long &N)
{
    //char hex_ch = '\0';
    unsigned long long *tmp = new unsigned long long[N.max_bit];
    int i = 0;
    int k = 0;
    unsigned long long elem = 0;
    unsigned long long tmp_el = 0;
    unsigned long long sum = 0;

    while (i < N.max_bit)
    {
        for (int j = i; j < i + 4; j++)
        {
            elem = N.bin[j];
            tmp_el = pow(2, j % 4) * elem;
            sum += tmp_el;
        }
        tmp[k] = sum;
        // cout<<tmp[k]<<endl;
        sum = 0;
        k++;
        i += 4;
    }

    unsigned long long hex_tmp = 0;
    for (int i = k - 1; i >= 0; i--)
    {
        hex_tmp = tmp[i];
        switch (hex_tmp)
        {
       case 0:
            N.hex_str = N.hex_str + '0';
            break;
       case 1:
            N.hex_str = N.hex_str + '1';
            break;
        case 2:
            N.hex_str = N.hex_str + '2';
            break;
        case 3:
            N.hex_str = N.hex_str + '3';
            break;
        case 4:
            N.hex_str = N.hex_str + '4';
            break;
        case 5:
            N.hex_str = N.hex_str + '5';
            break;
        case 6:
            N.hex_str = N.hex_str + '6';
            break;
        case 7:
            N.hex_str = N.hex_str + '7';
            break;
        case 8:
            N.hex_str = N.hex_str + '8';
            break;
        case 9:
            N.hex_str = N.hex_str + '9';
            break;
        case 10:
            N.hex_str = N.hex_str + 'A';
            break;
        case 11:
            N.hex_str = N.hex_str + 'B';
            break;
        case 12:
            N.hex_str = N.hex_str + 'C';
            break;
        case 13:
            N.hex_str = N.hex_str + 'D';
            break;
        case 14:
            N.hex_str = N.hex_str + 'E';
            break;
        case 15:
            N.hex_str = N.hex_str + 'F';
            break;
        
        default:
            break;
        }
    }
    
    cout << endl;    

}

unsigned long int MaxDeg(Long A)
{
    unsigned long int maxdeg = log2(A.coef[A.deg-1]) + 1 + (A.deg - 1) * 32;
    return maxdeg;
}



Long Add( Long A, Long B)
{
    Long C;
    C.deg = max(A.deg,B.deg);
    for (int i = 0; i < C.deg; i++){
        C.coef[i] = A.coef[i] ^ B.coef[i];
    }

    while (C.coef[C.deg-1]==0){
        C.deg--;
    }
    return C;
}


Long MultOneDig(Long A, long long b)
{
    Long C;
    long long base = 1;
    unsigned long long carry = 0;
    unsigned long long tmp = 0;
    base = (base << 32);
    for (int i = 0; i < A.deg; i++)
     {
        tmp = A.coef[i] * b + carry;
        C.coef[i] = tmp & (base - 1);
        C.deg++;
        carry = (tmp >> 32);
     }
    if(carry)
    {
        C.coef[C.deg] = carry;
        C.deg++;
    }
    return C;
}

Long ShiftBitsToHigh(Long N, int pos)
{
    for (int i = 0; i < pos; i++)
    {
        N = MultOneDig(N, 2);
    } 
    return N;
}
int Compare(Long A, Long B)
{
    if(A.deg > B.deg) 
        return 1;
    if(A.deg < B.deg)
        return -1;
    for (int i = A.deg - 1; i >= 0; i--)
    {
        if(A.coef[i] > B.coef[i])
        {
            return 1;
        }
        if (A.coef[i] < B.coef[i])
        {
            return -1;
        }
    }
    return 0;
}

Long ModGen(Long A)
{
    Long R = A;
    if(Gen.coef[Gen.deg - 1] != 0)
    {
        Long Tmp;
        Long Tmp1;
        int Gen_bit_len = MaxDeg(Gen);
        int R_bit_len;
        
        int cmp_R_Gen = Compare(R, Gen);
        int cmp_R_Tmp;
        while (cmp_R_Gen >= 0)
        {
            R_bit_len = MaxDeg(R);
            Tmp = ShiftBitsToHigh(Gen, R_bit_len - Gen_bit_len);
            cmp_R_Tmp = Compare(R, Tmp);
            while (cmp_R_Tmp == -1)
            {
                R_bit_len--;
                Tmp = ShiftBitsToHigh(Gen, R_bit_len - Gen_bit_len);
                cmp_R_Tmp = Compare(R, Tmp);
            }
            R = Add(R, Tmp);
            cmp_R_Gen = Compare(R, Gen);
        }
    }

    return R;
}

int Mask(Long N, int k)
{
    int mask = 1 << (k % 32);
    int  tmp_bit = k / 32;
    int answer = (N.coef[tmp_bit] & mask) >> (k % 32);
    return answer;
}


Long Mult(Long A, Long B)
{
    Long Res, TmpMult;
    int bitB=0;
    TmpMult = A;
    B.max_bit = MaxDeg(B);
    for (int i = 0; i < B.max_bit; i++)
    {
        bitB = Mask(B,i); // 
        if (bitB==1)    
        {
            Res = Add(Res, TmpMult);   
        }
        TmpMult = MultOneDig(TmpMult,2);        
    }

    while (Res.coef[Res.deg-1]==0)
    {
        Res.deg--;
    }
    Res = ModGen(Res);    
    return Res;
}


Long Square(Long A)
{
    Long N;
    for (int i = A.max_bit; i > 0; i--)
    {
        N.bin[2 * i] = A.bin[i];
    }
    N.max_bit = 2 * A.max_bit;
    bin_to_hex(N);
    Input(N);
    N = ModGen(N);
    return N;
}

Long Trace(Long A)
{
    Long Conj_tmp = A;
    Long Res;
    Long C;
    for (int i = 0; i < Gen.max_bit - 1; i++)
    {
        cout << i << endl;
        Conj_tmp = Square(Conj_tmp);
        Output_str(Conj_tmp);
        Write_bin(Conj_tmp);
        Input(Conj_tmp);
        Res = Add(Res, Conj_tmp);
    }
    return Res;  
}


Long Power(Long A, Long N)
{
    Write_bin(N);
    Write_bin(A);
    Long Res = A;
    for (int i = N.max_bit - 2; i >= 0; i--)
    {
        // cout << i << endl;
        Res = Mult(Res,Res);
        // Res = Square(Res);
        // Output_str(Res);
        // Write_bin(Res);
        // Input(Res);
        if (N.bin[i]==1)
        {
            Res = Mult(Res,A);
        }
    }
    
    return Res;

}


int main()
{

    unsigned int start_time = 0;
    unsigned int end_time = 0;
    double search_time = 0.0;
    Long A;
    Long B;
   
    Long Trace_res;
    Long Add_res;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
    Long Mult_res;
    Long Square_res;
    Long Power_res;
    Long Inv_res;
    Long Shift;
    
    Gen.bin[509]=1;
    Gen.bin[23]=1;
    Gen.bin[3]=1;
    Gen.bin[2]=1;
    Gen.bin[0]=1;
    Gen.max_bit=510;
    bin_to_hex(Gen);
    Input(Gen);
    A.hex_str = "180BC372B3C1166D5B0D879C0880B6B868B54860485834C9D2DDE25B6AF55DF3EFCEC2C908BFD5685C2C842BAD5798F2A7B2DC0A56028157CD2E8A288272E646";
    B.hex_str = "100F56108A35305C0ACC56DEACCA8A6A19AD474606B8617284AC51C3196AC1ECF6A7C1BF070362C8B7AAEC21A2CA1B057F2A63628ACFE1244A71C9EFEC96ABB5";
    
    Input(A);
    Input(B);

    
    
    


//VARIABLES:
    cout << "A:" << endl;
    Output(A);
    cout << "B:" << endl;
    Output(B);
// //Add
//     start_time = clock();
//     Add_res = Add(A,B);
//     end_time = clock();
//     cout << "\nAdd:" << endl;
//     Output(Add_res);
//     search_time = end_time - start_time;
//     cout << "Add time = " << search_time / 1000 <<endl;
// //Mult
//     start_time = clock();
//     Mult_res = Mult(A,B);
//     end_time = clock();
//     // Mult_res = ModGen(Mult_res);
//     cout << "\nMult:" << endl;
//     Output(Mult_res);
//     search_time = end_time - start_time;
//     cout << "Mult time = " << search_time / 1000 <<endl;

// //Power
//     Long M;
//     M.hex_str = "1F66C724D22924F9E89808F25C9986261A97EB379FFC1749580E33C4E8B696F213A4F2F2F65DC5394DBD1448B7C7866790DEE67490C925CD6CC02AB02D0B25D9";
//     Input(M);
//     start_time = clock();
//     Power_res = Power(A,M);
//     end_time = clock();
//     cout << "Power(A,M)" << endl;
//     Output(Power_res);
//     search_time = end_time - start_time;
//     cout << "Power time = " << search_time / 1000 <<endl;
    
// //Inverse
//     Long N; // N = 2^509 - 1
//     N.max_bit = 509;

//     for (int i = 1; i < N.max_bit; i++)
//     {
//         N.bin[i] = 1;
//     }
//     bin_to_hex(N);
//     Input(N);
    
//     start_time =clock();
//     Inv_res = Power(A,N);
//     end_time = clock();
//     search_time = (end_time - start_time) / 1000;
//     cout << "Inv(A^-1)" << endl;
//     Output(Inv_res);
//     cout << "Inv time = " << search_time << endl;
//     cout << endl;
//Square:
    Write_bin(A);
    start_time = clock();
    Square_res = Square(A);
    end_time = clock();
    cout << "\nSquare:" << endl;   
    Output(Square_res);
    search_time = end_time - start_time;
    cout << "Square time = " << search_time / 1000 <<endl;
// //Trace
//     Write_bin(A);
//     Trace_res = Trace(A);
//     bin_to_hex(Trace_res);
    return 0;   
}