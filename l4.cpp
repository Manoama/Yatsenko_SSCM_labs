#include<iostream>
#include<cmath>
#include<string>
#include <typeinfo>
#include <ctime>

using namespace std;

struct Long
{
    string hex_str;
    int deg=0;
    long long coef[33] = { 0 };
    int max_bit=0;
    int bin[2999] = { 0 };
};


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

Long Add(Long A, Long B){
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


//Цикличекий сдвиг
Long Square(Long A)
{
    Long C;
    for (int i = 0; i < A.max_bit - 1; i++){
        C.bin[i] = A.bin[i + 1];
    }
    C.bin[172] = A.bin[0];
    if (C.bin[172]) // ==1   
        C.max_bit = 173;
    else
        C.max_bit = A.max_bit - 1;
    return C;
}

int Trace(Long A)
{
    int trace;
    trace = 0;
    for (int i = 0; i <= (173 / 3); i++){
        trace = trace ^(A.bin[3 * i] ^ A.bin[3 * i + 1] ^ A.bin[3 * i + 2]);
    }
    return trace;    
}

// Матрица Л

// 2 ^ n mod 347(2*173 + 1)
int power_2(int n)
{
    int sum = 1;
    int p = 2 * 173 + 1;
    for (int i = 0; i < n; i++)
    {
        sum *= 2;
        sum = sum % p;
    }
    
    return sum;

}


// used in Mult()
void ShiftBin(Long &A)
{
    Long Tmp = A;
    for (int i = 1; i < 173; i++)
    {
        Tmp.bin[i] = A.bin[i-1];
    }
    Tmp.bin[0] = A.bin[172];
    for (int i = 0; i < 173; i++)
    {
        A.bin[i] = Tmp.bin[i];
    }
}

// multiplicative Matrix M0
void Matrix(int M[173][173])
{
    int elem_i = 0 ;
    int elem_j = 0;
    int a,b;
    int p = 2 * 173 + 1;
    for (int i = 0; i < 173; i++){
        for (int j = 0; j < 173; j++){
            elem_i = power_2(i);
            elem_j = power_2(j);
            a = (elem_i + elem_j) % p;
            b = (elem_i - elem_j) % p;
    
            if(a == 1 || b == 1 || a == -1 || b == -1 || a == -346 || b == -346 || a == 346 || b == 346)
                M[i][j] = 1;
        }   
    }
}

//Mult
Long Mult(Long A, Long B, int M[173][173])
{
    Long Res;
    int column;
    int b = 0;
    int Tmp[173] = { 0 };
    int sum = 0;
    int mult = 0;
    int bit = 0;
    int res = 0;
    for (int l = 0; l < 173; l++){
        for (int j = 0; j < 173; j++){
            for (int i = 0; i < 173; i++){
                bit = M[i][j];
                mult = bit * A.bin[172 - i];
                sum += mult;
            }
            column = sum % 2;
            Tmp[j] = column;
            bit = 0;
            sum = 0;
            mult = 0;
        }
        res = 0;
        for (int k = 0; k < 173; k++){
            res = (res ^ (Tmp[k]*B.bin[172 - k]));
        }      
        ShiftBin(A);
        ShiftBin(B);
        Res.bin[172-l] = res;
    }
    Res.max_bit = 173;
    return Res;
}

Long Power(Long A, Long N, int M[173][173])
{
    Long Res;
    Res = A;
    for (int i = N.max_bit - 2; i >= 0; i--){
        Res = Square(Res);
        if (N.bin[i])
            Res = Mult(Res,A,M);
    }
    return Res;
}

// m - 1 = 172
Long Inv(Long A, int M[173][173])
{
    
    Long B_curr = A;
    Long B_prev = A;
    int k = 1;
    // 172 = 10101100
    int m[8];
    m[0] = 1;
    m[1] = 0;
    m[2] = 1;
    m[3] = 0;
    m[4] = 1;
    m[5] = 1;
    m[6] = 0;
    m[7] = 0; 
    for (int i = 1; i < 8; i++){
        for (int j = 0; j < k; j++){
            B_curr = Square(B_curr);
        }
        B_curr = Mult(B_curr, B_prev, M);
        k *= 2;
        if (m[i]==1){
            B_curr = Square(B_curr);
            B_curr = Mult(B_curr, A, M);
            k++;
        }
        B_prev = B_curr;        
    }
    B_curr = Square(B_curr);
    return B_curr;
}

int main(void)
{
    Long A,B,N;
    Long Add_res;
    Long AddBin_res;
    Long Square_res;
    Long Mult_res;
    Long Power_res;
    Long Inv_res;
    A.hex_str="198033D9EBC41A593612C054B6336C856952DB545DC9";
    B.hex_str="1320C5786E1906CCD14EC313BD9D0330870FD35A6884";
    N.hex_str="1AAD6013A54EFD7E984D8F52502E849CEE90D74E1F31";
    Input(A);
    Input(B);
    Write_bin(A);
    Write_bin(B);

    int M[173][173];
    Matrix(M);

    unsigned int start_c = 0;
    unsigned int end_c = 0;
    double search_time = 0.0;

// //Addition
//     start_c = clock();
//     Add_res = Add(A,B);
//     end_c = clock();
//     cout << "\nAdd:" << endl;
//     Output(Add_res);
//     search_time = end_c - start_c;
//     cout << "Add time = " << search_time / 1000 << endl;
// // Square
//     Write_bin(A);
//     start_c = clock();
//     Square_res = Square(A);
//     end_c = clock();
//     bin_to_hex(Square_res);
//     cout << "\nSquare_res = " <<Square_res.hex_str << endl;   
//     search_time = end_c - start_c;
//     cout << "Add time = " << search_time / 1000 << endl; 
// //Mult
//     start_c = clock();
//     Mult_res = Mult(A,B,M);
//     end_c = clock();
//     bin_to_hex(Mult_res);
//     cout << "\nMult_res = " << Mult_res.hex_str << endl;
//     search_time = end_c - start_c;
//     cout << "Add time = " << search_time / 1000 << endl;
// // Trace(A)
//     Write_bin(A);
//     start_c = clock();
//     cout << "\nTrace(A) = " << Trace(A) << endl;
//     end_c = clock();
//     search_time = end_c - start_c;
//     cout << "Add time = " << search_time / 1000 << endl;
// // Trace(B)
//     Write_bin(B);
//     start_c = clock();
//     cout << "\nTrace(B) = " << Trace(B) << endl;
//     end_c = clock();
//     search_time = end_c - start_c;
//     cout << "Add time = " << search_time / 1000 << endl;

// //Power
//     Write_bin(A);
//     Write_bin(N);
//     start_c = clock();
//     Power_res = Power(A,N,M);
//     end_c = clock();
//     bin_to_hex(Power_res);
//     cout << "\nPower_res = " << Power_res.hex_str << endl;
//     search_time = end_c - start_c;
//     cout << "Add time = " << search_time / 1000 << endl;
// //Inv
//     Write_bin(A);
//     start_c = clock();
//     Inv_res = Inv(A,M);
//     end_c = clock();
//     bin_to_hex(Inv_res);
//     cout << "\nInv_res = " << Inv_res.hex_str << endl;
//     search_time = end_c - start_c;
//     cout << "Add time = " << search_time / 1000 << endl;
    






















    return 0;
}