#include "kbhash2.h"
#include <time.h>
using namespace std;


string kbhash::int2kmer(long long x, int k)
{
    string s;

    for(int i = 0; i < k; i++)
    {
        int y = x & 3;
        x = (x>>2);

        if(y == 0)
            s = 'A' + s;
        if(y == 1)
            s = 'C' + s;
        if(y == 2)
            s = 'G' + s;
        if(y == 3)
            s = 'T' + s;
    }

    return s;
}

long long kbhash::kmer2int(string s)
{
    long long ret = 0;
    int k = s.length();

    for(int i = 0; i < k; i++)
    {
        if(s[i] == 'A')
            ret += 0;
        if(s[i] == 'C')
            ret += 1;
        if(s[i] == 'G')
            ret += 2;
        if(s[i] == 'T')
            ret += 3;

        if(i < k - 1)
            ret <<= 2;
    }

    return ret;
}


kbhash::kbhash(int k, int b, int num, int m)
{
    int k1 = k - b;
    klen = k;
    blen = b;
    vlen = num;
    mlen = m;
    clock_t t0, t1, t2;
    t0 = clock();

    // ofstream outfile;
    // outfile.open("./bmers");

    for(int i = 0; i < num; i++)
    {
        long long x = 0;
        for(int j = 0; j < 2*k1; j++)
            x = x + (((long long)(rd() % 2)) << j);

        ranks.push_back(x);
        // string km = int2kmer(x, k1);
        // outfile<<x<<endl;

        //bestb.push_back(km);
    }

    t1 = clock();
    dict['A'] = 0;
    dict['C'] = 1;
    dict['G'] = 2;
    dict['T'] = 3;


    build_hash_table(m);
    t2 = clock();

    printf("%.2lf\t%.2lf\n", double(t1 - t0)/CLOCKS_PER_SEC, double(t2 - t1)/CLOCKS_PER_SEC);
}

kbhash::kbhash(int k, int b, vector<long long> bmers, int m)
{
    int k1 = k - b;
    klen = k;
    blen = b;
    vlen = bmers.size();
    mlen = m;

    for(int i = 0; i < vlen; i++)
    {
        long long x = bmers[i];

        ranks.push_back(x);
        // string km = int2kmer(x, k1);
        // bestb.push_back(km);
    }

    dict['A'] = 0;
    dict['C'] = 1;
    dict['G'] = 2;
    dict['T'] = 3;


    build_hash_table(m);
}

void kbhash::build_hash_table(int m)
{
    int num = (1<<(2 * m));
    // int paths = (1<<(2 * numt - 2));

    for(int i = 0; i < num; i++)
        heads.push_back(0);

    int k1 = klen - blen;

    for(int i = 0; i < vlen; i++)
    {
        long long id = ranks[i];

        int v = 0;
        // int r = 0;
        // int l = 0;
        // int sumc[4] = {0};
        // for(int j = 0; j < klen - blen; j++)
        //     sumc[dict[s[j]]]++;

        // int key = 0, maxi = sumc[0];
        // for(int j = 1; j < 4; j++)
        //     if(sumc[j] > maxi)
        //     {
        //         key = j;
        //         maxi = sumc[j];
        //     }

        v = (id >> ((k1-m)*2));
        // for(int j = m; j < m + numt; j++)
        //     r = (r<<2) + dict[s[j]];
        

        hashtable* tmp = new hashtable;
        tmp->rank = i;
        tmp->pos = 0;
        tmp->next = heads[v];
        // tmp->lnum = -1;
        // tmp->rnum = r;

        // for(int k = 0; k < 4; k++)
        //     tmp->sum[k] = sumc[k];
        // tmp->key = key;
        // tmp->num = maxi;

        heads[v] = tmp;
        
        for(int j = m; j < k1; j++)
        {
            v = v % (num>>2);
            v = (v<<2) + ((id>>((k1 - j - 1) * 2)) & 3);


            // l = l % paths;
            // l = (l<<2) + dict[s[j-m]];


            hashtable* tmp = new hashtable;
            tmp->rank = i;
            tmp->pos = j - m + 1;
            tmp->next = heads[v];

            // if(j + numt < klen - blen)
            // {
            //     r = r % paths;
            //     r = (r<<2) + dict[s[j+numt]];
            //     tmp->rnum = r;
            // }
            // else
            //     tmp->rnum = -1;

            // if(j - m >= numt - 1)
            //     tmp->lnum = l;
            // else
            //     tmp->lnum = -1;
            // for(int k = 0; k < 4; k++)
            //     tmp->sum[k] = sumc[k];
            // tmp->key = key;
            // tmp->num = maxi;

            heads[v] = tmp;
        }
    }
}

readinfo* kbhash::get_Bneighbor(string s)
{       
    int len = s.length();
    int k1 = klen - blen;
    readinfo* r = new readinfo(len, klen, k1);

    int next[4] = {-1, -1, -1, -1};

    for(int i = 0; i < 4; i++)
        edges[len-1][i] = -1;

    for(int i = len - 1; i > 0; i--)
    {
        int c = dict[s[i]];
        next[c] = i;

        for(int j = 0; j < 4; j++)
            edges[i-1][j] = next[j];
    }

    int num = len - klen;
    next[dict[s[0]]] = 0;

    if(len < klen)
        return r;

    deque<int> indices;

    for(int i = 0; i < vlen; i++)
    {
        long long x = ranks[i];
        int firstc = (x>>(2 * (k1 - 1))) & 3;
        int st = next[firstc];
        

            
        while(st <= num && st != -1)
        {
            int state = st;

            indices.clear();
            indices.push_back(st);

            for(int j = k1 - 2; j >= 0; j--)
            {
                int y = (x>>(2 * j)) & 3;
                // cout<<state<<" "<<y<<endl;
                // cout<<edges[state][y]<<endl;
                state = edges[state][y];
                indices.push_back(state);

                if(state - st > klen || state == -1)
                    break;
            }

            if(state != -1 && state - st < klen)
            {
                bnbr tmp;
                tmp.st = st;
                tmp.en = state;
                tmp.id = i;
                for(int index: indices)
                    tmp.indices.push_back(index);

                //cout<<st<<" "<<state<<" "<<i<<endl;
                if(r->slist[st][state-st+1-k1].id == -1 || r->slist[st][state-st+1-k1].id > i)
                {
                    // cout<<st<<" "<<state<<" "<<state-st+1-k1<<" "<<i<<endl;
                    r->slist[st][state-st+1-k1] = tmp;
                }
            }

            st = edges[st][firstc];
        }
    }

    r->getminbnlist();
    return r;
}

void kbhash::concat(readinfo* r, string s)
{
    r->getminbnlist();

    vector<bnbr> newbn;
    for(auto b: r->minbnlist)
    {
        for(int i = b.vst; i <= b.ven - klen + 1; i++)
            if(i >= flen)
            {
                bnbr tmp;
                tmp.st = i - flen;
                tmp.en = i + klen - 1;
                tmp.id = b.id + vlen * kmer2int(s.substr(tmp.st, flen));

                //cout<<tmp.st<<" "<<tmp.en<<" "<<tmp.id<<endl;
                newbn.push_back(tmp);
            }
    }

    r->minbnlist = newbn;
}
readinfo* kbhash::get_Bneighbor_subseq(string s)
{       
    clock_t t0, t1, t2;
    t0 = clock();

    int len = s.length();
    int k1 = klen - blen;
    readinfo* r = new readinfo(len, klen, k1);

    int next[4] = {-1, -1, -1, -1};
    int last[4] = {-1, -1, -1, -1};



    for(int i = 0; i < 4; i++)
    {
        edges[len-1][i] = -1;
        edgesv[0][i] = -1;
        //sum[0][i] = 0;
    }

   
    // for(int i = 1; i <= len; i++)
    // {        
    //     for(int j = 0; j < 4; j++)
    //         sum[i][j] = sum[i-1][j];
    //     sum[i][dict[s[i-1]]]++;
    // }

    for(int i = len - 1; i > 0; i--)
    {
        int c = dict[s[i]];
        next[c] = i;

        for(int j = 0; j < 4; j++)
            edges[i-1][j] = next[j];
    }

    for(int i = 0; i < len - 1; i++)
    {
        int c = dict[s[i]];
        last[c] = i;

        for(int j = 0; j < 4; j++)
            edgesv[i+1][j] = last[j];
    }

    int num = len - klen + 1;
    int mnum = len - mlen + 1;

    if(num <= 0)
        return r;

    int indices[50];
    

    int x = 0;
    for(int i = 0; i < mlen - 1; i++)
        x = (x<<2) + dict[s[i]];

    int lnum = 0;
    int l2num = 0;
    int rnum = 0;

    for(int i = 0; i < mnum; i++)
    {
        x = x % (1<<(2*mlen-2));
        x = (x<<2) + dict[s[i + mlen - 1]];

        for(hashtable* j = heads[x]; j != NULL; j = j->next)
        {
            int id = j->rank;
            int pos = j->pos;

            int st = i, en = i + mlen - 1;
            int tst = i, ten = i + mlen - 1;

            lnum++;
            if(st < pos || en + k1 - pos - mlen >= len)
                continue;

            // int lm = max(0, st - pos - blen);
            // int rm = min(len, en + k1 - pos - mlen + blen + 1);

            // //cout<<lm<<" "<<rm<<endl;

            //     // cout<<bestb[id]<<endl;
            //     // cout<<j->sum[k]<<" "<<sum[rm][k]<<" "<<sum[lm][k]<<endl;
            // if(j->num > sum[rm][j->key] - sum[lm][j->key])
            //     continue;

            long long nbr = ranks[id]; 
            l2num++;
            //cout<<j<<" "<<h<<" "<<id<<endl;
            // cout<<i<<" "<<pos<<" "<<bestb[id]<<endl;
            // cout<<j->lnum<<" "<<j->rnum<<endl;

            int forw = pos - 1, back = pos + mlen, y;

            // cout<<st<<" "<<en<<endl;
            while(forw >= 0 || back < k1)
            {   
                if(forw >= 0)
                {
                    y = (nbr>>(2 * (k1 - forw - 1))) & 3;
                    st = edgesv[st][y];

                    if(st == -1)
                        break;

                    if(en + forw + (k1 - back) - st >= klen)
                    {
                        st = -1;
                        break;
                    }

                    forw--;
                }

                if(back < k1)
                {
                    y = (nbr>>(2 * (k1 - back - 1))) & 3;
                    en = edges[en][y];

                    if(en == -1)
                        break;     

                    if(en + forw + (k1 - back) - st >= klen)
                    {
                        en = -1;
                        break;
                    }
        
                    back++;
                }
            }

            // cout<<i<<" "<<id<<" "<<st<<" "<<en<<endl;
            // cout<<bestb[id]<<" "<<pos<<endl;
            
            if(st == -1 || en == -1 || en - st + 1 > klen)
                continue;

            bnbr tmp;
            tmp.st = st;
            tmp.en = en;
            tmp.est = i;
            tmp.een = i + mlen - 1; 
            rnum++;

            tmp.id = id;

            int nst = st;

            for(int k = 1; k < pos; k++)
            {
                tmp.indices.push_back(nst);
                nst = edges[nst][(nbr>>(2 * (k1 - k - 1))) & 3];
            }
            if(st < tst)
                tmp.indices.push_back(nst);
            
            for(int k = tst; k <= ten; k++)
                tmp.indices.push_back(k);

            nst = ten;
            for(int k = pos + mlen; k < k1; k++)
            {
                nst = edges[nst][(nbr>>(2 * (k1 - k - 1))) & 3];
                tmp.indices.push_back(nst);
            }


//            cout<<pos<<endl;
            // cout<<bestb[id]<<endl;
            // for(auto ind:tmp.indices)
            //     cout<<ind<<" ";
            // cout<<endl;

            if(r->slist[st][en-st+1-k1].id == -1 || r->slist[st][en-st+1-k1].id > id)
                r->slist[st][en-st+1-k1] = tmp;
        }
    }

    t1 = clock();
    r->getminbnlist();
    t2 = clock();
    //concat(r, s);
    //printf("%.2lf\t%.2lf\n", double(t1 - t0)/CLOCKS_PER_SEC, double(t2 - t1)/CLOCKS_PER_SEC);
    ///printf("%.2lf\t%.2lf\n", double(lnum)/mnum, double(l2num)/mnum, double(rnum)/mnum);
    return r;
}

// readinfo* kbhash::get_Bneighbor_subseq2(string s)
// {       
//     clock_t t0, t1, t2;
//     t0 = clock();

//     int len = s.length();
//     int k1 = klen - blen;
//     readinfo* r = new readinfo(len, klen, k1);

//     int next[4] = {-1, -1, -1, -1};
//     int last[4] = {-1, -1, -1, -1};



//     for(int i = 0; i < 4; i++)
//     {
//         edges[len-1][i] = -1;
//         edgesv[0][i] = -1;
//         //sum[0][i] = 0;
//     }

   
//     // for(int i = 1; i <= len; i++)
//     // {        
//     //     for(int j = 0; j < 4; j++)
//     //         sum[i][j] = sum[i-1][j];
//     //     sum[i][dict[s[i-1]]]++;
//     // }
//     for(int i = 0; i < (1<<(numt * 2)); i++)
//         pathl[i] = pathr[i] = -1;

//     for(int i = len - 1; i > 0; i--)
//     {
//         int c = dict[s[i]];
//         next[c] = i;

//         for(int j = 0; j < 4; j++)
//             edges[i-1][j] = next[j];
//     }

//     for(int i = 0; i < len - 1; i++)
//     {
//         int c = dict[s[i]];
//         last[c] = i;

//         for(int j = 0; j < 4; j++)
//             edgesv[i+1][j] = last[j];
//     }

//     int num = len - klen + 1;
//     int mnum = len - mlen + 1;

//     if(num <= 0)
//         return r;

//     int indices[50];
    

//     int x = 0;
//     for(int i = 0; i < mlen - 1; i++)
//         x = (x<<2) + dict[s[i]];

//     int lnum = 0;
//     int l2num = 0;
//     int rnum = 0;

//     for(int i = 0; i < mnum; i++)
//     {
//         x = x % (1<<(2*mlen-2));
//         x = (x<<2) + dict[s[i + mlen - 1]];

//         int z;
//         if(i > 0)
//         {
//             z = dict[s[i-1]];
//             for(int k = z ; k < (1<<(numt * 2)); k += 4)
//                 pathl[k] = -1;
//         }

//         z = dict[s[i+mlen-1]];
//         for(int k = z * (1<<(numt * 2 - 2)); k < (z+1) * (1<<(numt * 2 - 2)); k++)
//             pathr[k] = -1;

//         //cout<<i<<" "<<(z * (1<<(numt * 2 - 2)))<<" "<<(z+1) * (1<<(numt * 2 - 2))<<endl;

//         for(hashtable* j = heads[x]; j != NULL; j = j->next)
//         {
//             int id = j->rank;
//             int pos = j->pos;

//             int st = i, en = i + mlen - 1;
//             int tst = i, ten = i + mlen - 1;

//             lnum++;
//             if(st < pos || en + k1 - pos - mlen >= len)
//                 continue;

//             // int lm = max(0, st - pos - blen);
//             // int rm = min(len, en + k1 - pos - mlen + blen + 1);

//             // //cout<<lm<<" "<<rm<<endl;

//             //     // cout<<bestb[id]<<endl;
//             //     // cout<<j->sum[k]<<" "<<sum[rm][k]<<" "<<sum[lm][k]<<endl;
//             // if(j->num > sum[rm][j->key] - sum[lm][j->key])
//             //     continue;

//             long long nbr = ranks[id]; 
//             l2num++;
//             //cout<<j<<" "<<h<<" "<<id<<endl;
//             // cout<<i<<" "<<pos<<" "<<id<<" "<<bestb[id]<<endl;
//             // cout<<j->lnum<<" "<<j->rnum<<endl;

//             int forw = pos - 1, back = pos + mlen, y;
//             if(j->lnum > 0)
//                 if(pathl[j->lnum] != -1)
//                 {
//                     // cout<<pathl[j->lnum]<<endl;
//                     forw -= numt;
//                     st = pathl[j->lnum];
//                 }

//             if(j->rnum > 0)
//                 if(pathr[j->rnum] != -1)
//                 {
//                     // cout<<pathr[j->rnum]<<endl;
//                     back += numt;
//                     en = pathr[j->rnum];
//                 }

//             // cout<<st<<" "<<en<<endl;
//             while(forw >= 0 || back < k1)
//             {   
//                 if(forw >= 0)
//                 {
//                     y = (nbr>>(2 * (k1 - forw - 1))) & 3;
//                     st = edgesv[st][y];
//                     // cout<<"forw: "<<y<<" "<<st<<endl;
//                     if(en + forw + (k1 - back) - st >= klen)
//                     {
//                         st = -1;
//                         break;
//                     }

//                     if(st == -1)
//                         break;

//                     if(forw == pos - numt)
//                         pathl[j->lnum] = st;

//                     forw--;
//                 }

//                 if(back < k1)
//                 {
//                     y = (nbr>>(2 * (k1 - back - 1))) & 3;
//                     en = edges[en][y];
//                     // cout<<"back: "<<y<<" "<<en<<endl;

//                     if(en + forw + 1 + (k1 - back - 1) - st >= klen)
//                     {
//                         en = -1;
//                         break;
//                     }


//                     if(en == -1)
//                         break;               

//                     if(back == pos + mlen + numt - 1)
//                         pathr[j->rnum] = en;

//                     back++;
//                 }
//             }
//             // cout<<forw<<" "<<back<<endl;
//             // cout<<st<<" "<<en<<endl;
//             if(st == -1 || en == -1 || en - st + 1 > klen)
//                 continue;

//             bnbr tmp;
//             tmp.st = st;
//             tmp.en = en;
//             tmp.est = i;
//             tmp.een = i + mlen - 1; 
//             rnum++;

//             tmp.id = id;

//             int nst = st;

//             for(int k = 1; k < pos; k++)
//             {
//                 tmp.indices.push_back(nst);
//                 nst = edges[nst][(nbr>>(2 * (k1 - k - 1))) & 3];
//             }
//             if(st < tst)
//                 tmp.indices.push_back(nst);
            
//             for(int k = tst; k <= ten; k++)
//                 tmp.indices.push_back(k);

//             nst = ten;
//             for(int k = pos + mlen; k < k1; k++)
//             {
//                 nst = edges[nst][(nbr>>(2 * (k1 - k - 1))) & 3];
//                 tmp.indices.push_back(nst);
//             }


// //            cout<<pos<<endl;
//             // cout<<id<<" "<<bestb[id]<<endl;
//             // for(auto ind:tmp.indices)
//             //     cout<<ind<<" ";
//             // cout<<endl;

//             if(r->slist[st][en-st+1-k1].id == -1 || r->slist[st][en-st+1-k1].id > id)
//                 r->slist[st][en-st+1-k1] = tmp;
//         }
//     }

//     t1 = clock();
//     r->getminbnlist();
//     t2 = clock();
//     //concat(r, s);
//     //printf("%.2lf\t%.2lf\n", double(t1 - t0)/CLOCKS_PER_SEC, double(t2 - t1)/CLOCKS_PER_SEC);
//     ///printf("%.2lf\t%.2lf\n", double(lnum)/mnum, double(l2num)/mnum, double(rnum)/mnum);
//     return r;
// }

vector<int> kbhash::suffix_function(string s) 
{
    int n = (int)s.length();
    vector<int> pi(n, 0);

    for (int i = 1; i < n; i++)
    {
        int j = pi[i-1];
        while(j && s[n-i-1] != s[n-j-1])
            j = pi[j-1];

        if(s[n-i-1] == s[n-j-1])
            j++;

        pi[i] = j;
    }

    return pi;
}

vector<int> kbhash::prefix_function(string s) 
{
    int n = (int)s.length();
    vector<int> pi(n, 0);

    for (int i = 1; i < n; i++)
    {
        int j = pi[i-1];
        while(j && s[i] != s[j])
            j = pi[j-1];

        if(s[i] == s[j])
            j++;

        pi[i] = j;
    }

    return pi;
}


segment_gap kbhash::prematch(string t, string p)
{
    vector<int> ppi = prefix_function(p);
    vector<int> spi = suffix_function(p);

    int lent = t.length();
    int lenp = p.length();

    int i = 0;
    int k = 0;
    int x[200];
    int y[200];
    int pos[200];

    for(int j = 0; j < lent; j++)
    {
        x[j] = -1;
        y[j] = lenp;
        pos[j] = j;
    }

    while(i < lent)
    {
        if(p[k] == t[i])
        {
            k++;
            i++;

            if(k == lenp)
            {
                x[i-1] = k-1;
                k = ppi[k-1];
            }
        }
        else
            if(k == 0)
                i++;
            else
                k = ppi[k - 1];

        if(i != 0 && k-1 > x[i-1])
            x[i-1] = k-1;
    }   

    i = lent - 1;
    k = lenp - 1;

    while(i >= 0)
    {
        if(p[k] == t[i])
        {
            k--;
            i--;

            if(k == -1)
            {
                y[i+1] = k+1;
                k = lenp - spi[lenp - 1] - 1;
            }
        }
        else
            if(k == lenp - 1)
                i--;
            else
                k = lenp - spi[lenp - k - 2] - 1;

        if(i != lent - 1 && k+1 < y[i+1])
            y[i+1] = k+1;
    }


    // for(i = 1; i < lent; i++)
    //  x[i] = max(x[i], x[i-1]);

    for(i = lent - 2; i >= 0; i--)
        if(y[i] > y[i+1])
        {
            y[i] = y[i+1];
            pos[i] = pos[i+1];
        }

    // cout<<t<<endl;
    // cout<<p<<endl;
    // for(int i = 0; i < lent - 1; i++)
    //  cout<<i<<" "<<x[i]<<endl;
    for(i = lent - 2; i >= 0; i--)
        if(x[i] >= y[i+1] - 1)
        {
            // cout<<t<<endl;
            // cout<<p<<endl;
            // cout<<i<<" "<<x[i]<<" "<<i+1<<" "<<y[i+1]<<" "<<pos[i+1]<<endl;
            //cout<<i-x[i]<<" "<<i-x[i]+y[i+1]-1<<" "<<pos[i+1]<<" "<<pos[i+1] + lenp - y[i+1] - 1<<endl;
            //return i - x[i];

            if(x[i] == -1)
                return (segment_gap){pos[i+1], pos[i+1] + lenp - y[i+1] - 1, -1, -1};

            if(y[i+1] == lenp) 
                return (segment_gap){i - x[i], i, -1, -1};
            
            return (segment_gap){i - x[i], i-x[i]+y[i+1]-1, pos[i+1], pos[i+1] + lenp - y[i+1] - 1};
        }

    //return -1;
    return (segment_gap){-1, -1, -1, -1};
}

segment_gap kbhash::sufmatch(string t, string p)
{
    vector<int> ppi = prefix_function(p);
    vector<int> spi = suffix_function(p);

    int lent = t.length();
    int lenp = p.length();

    int i = 0;
    int k = 0;
    int x[200];
    int y[200];
    int pos[200];

    for(int j = 0; j < lent; j++)
    {
        x[j] = -1;
        y[j] = lenp;
        pos[j] = j;
    }

    while(i < lent)
    {
        if(p[k] == t[i])
        {
            k++;
            i++;

            if(k == lenp)
            {
                x[i-1] = k-1;
                k = ppi[k-1];
            }
        }
        else
            if(k == 0)
                i++;
            else
                k = ppi[k - 1];

        if(i != 0 && k-1 > x[i-1])
            x[i-1] = k-1;
    }   

    i = lent - 1;
    k = lenp - 1;

    while(i >= 0)
    {
        if(p[k] == t[i])
        {
            k--;
            i--;

            if(k == -1)
            {
                y[i+1] = k+1;
                k = lenp - spi[lenp - 1] - 1;
            }
        }
        else
            if(k == lenp - 1)
                i--;
            else
                k = lenp - spi[lenp - k - 2] - 1;

        if(i != lent - 1 && k+1 < y[i+1])
            y[i+1] = k+1;
    }

    for(i = 1; i < lent; i++)
        if(x[i] < x[i-1])
        {
            x[i] = x[i-1];
            pos[i] = pos[i-1];
        }

    // for(i = lent - 2; i >= 0; i--)
    //  y[i] = min(y[i], y[i+1]);

    for(i = 0; i < lent - 1; i++)
        if(x[i] >= y[i+1] - 1)
        {           
            // cout<<t<<endl;
            // cout<<p<<endl;
            // cout<<pos[i] - x[i]<<" "<<pos[i]<<" "<<i+1<<" "<<i + lenp - x[i] - 1<<endl;
            //return i + lenp - y[i+1];

            if(x[i] == -1)
                return (segment_gap){i+1, i + lenp - y[i+1], -1, -1};

            if(y[i+1] == lenp) 
                return (segment_gap){pos[i] - x[i], pos[i], -1, -1};
            
            return (segment_gap){pos[i] - x[i], pos[i], i+1, i + lenp - y[i+1]};
        }

    //return -1;
    return (segment_gap){-1, -1, -1, -1};

    //return -1;
}

readinfo* kbhash::get_Bneighbor_substr4gap(string s)
{       
    int len = s.length();
    int k1 = klen - blen;

    readinfo* r = new readinfo(len, klen, k1);



    int next[4] = {-1, -1, -1, -1};
    int last[4] = {-1, -1, -1, -1};

    for(int i = 0; i < 4; i++)
    {
        edges[len-1][i] = -1;
        edgesv[0][i] = -1;
    }

    for(int i = len - 1; i > 0; i--)
    {
        int c = dict[s[i]];
        next[c] = i;

        for(int j = 0; j < 4; j++)
            edges[i-1][j] = next[j];
    }

    for(int i = 0; i < len - 1; i++)
    {
        int c = dict[s[i]];
        last[c] = i;

        for(int j = 0; j < 4; j++)
            edgesv[i+1][j] = last[j];
    }

    int num = len - klen + 1;
    int mnum = len - mlen + 1;

    if(num <= 0)
        return r;
    vector<int> indices;
    int x = 0;
    for(int i = 0; i < mlen - 1; i++)
        x = (x<<2) + dict[s[i]];

    for(int i = 0; i < mnum; i++)
    {
        x = x % (1<<(2*mlen-2));
        x = (x<<2) + dict[s[i + mlen - 1]];

        for(hashtable* j = heads[x]; j != NULL; j = j->next)
        {
            int id = j->rank;
            int pos = j->pos;
            int st = i, en = i + mlen - 1;
            int nbr = ranks[id]; 
            string c = bestb[id];

            int stext = pos - 1;
            int enext = pos + mlen;

            int exmatchst = i;
            int exmatchen = i + mlen - 1;
            
            indices.clear();
            for(int k = st; k <= en; k++)
                indices.push_back(k);

            if(pos != 0)
            {
                for( ;stext >= 0; stext--)
                {
                    int y = (nbr>>(2 * (k1 - stext - 1))) & 3;
                    if(st != edgesv[st][y] + 1) 
                        break;

                    exmatchst--;
                    st = edgesv[st][y];
                    indices.push_back(st);

                    if(st == -1)
                        break;
                }

                if(stext >= 0)
                {   
                    st--;

                    if(st >= 0)
                    {
                        string ptr = c.substr(0, stext + 1);


                        int tst = max(0, st - (blen + stext));
                        int ten = st;

                        segment_gap mpos = prematch(s.substr(tst, ten - tst + 1), ptr);

                        if(mpos.st1 == -1)
                            st = -1;
                        else
                            st = tst + mpos.st1;                        

                        if(mpos.st1 != -1)
                            for(int k = mpos.st1; k <= mpos.en1; k++)
                                indices.push_back(tst + k);
                        if(mpos.st2 != -1)
                            for(int k = mpos.st2; k <= mpos.en2; k++)
                                indices.push_back(tst + k);
                    }
                }
            }

            if(pos + mlen != k1)
            {
                for(; enext < k1; enext++)
                {
                    int y = (nbr>>(2 * (k1 - enext - 1))) & 3;
                    if(en != edges[en][y] - 1) 
                        break;

                    en = edges[en][y];
                    exmatchen++;
                    indices.push_back(en);

                    if(en == -1)
                        break;
                }

                if(enext < k1 && en != -1)
                {
                    en++;
                    if(en < len)
                    {
                        string ptr = c.substr(enext);

                        int ten = min(len - 1, en + (blen + k1 - enext) - 1);
                        int tst = en;

                        segment_gap mpos = sufmatch(s.substr(tst, ten - tst + 1), ptr);
                        
                        if(mpos.en1 == -1)
                            en = -1;
                        else
                        {   

                            if(mpos.st2 == -1)
                                en = tst + mpos.en1;
                            else
                                en = tst + mpos.en2;

                        }                       

                        if(mpos.st1 != -1)
                            for(int k = mpos.st1; k <= mpos.en1; k++)
                                indices.push_back(tst + k);
                        if(mpos.st2 != -1)
                            for(int k = mpos.st2; k <= mpos.en2; k++)
                                indices.push_back(tst + k);

                    }
                    else
                        en = -1;
                }
            }

            if(st < 0 || en < 0 || en - st + 1 > klen)
                continue;

            // cout<<c<<endl;
            // cout<<s.substr(st, en - st + 1)<<endl;
            // // if(pos != 0)
            // //   cout<<c.substr(0, stext + 1)<<" ";
            // // cout<<c.substr(stext + 1, enext - stext - 1)<<" ";
            // // if(pos + mlen != k1)
            // //   cout<<c.substr(enext);
            // // cout<<endl;
            // for(int i = 0; i < 5; i++)
            //  printf("[%d ,%d]\t", segs[i].first, segs[i].second);
            // printf("\n\n");

            bnbr tmp;
            tmp.st = st;
            tmp.en = en;
            tmp.est = i;
            tmp.een = i + mlen - 1;

            for(int index:indices)
                tmp.indices.push_back(index);

            tmp.id = id;

            if(r->slist[st][en-st+1-k1].id == -1 || r->slist[st][en-st+1-k1].id > id)
                r->slist[st][en-st+1-k1] = tmp;            
        }
    }

    r->getminbnlist();
    //concat(r, s);
    return r;
}


readinfo* kbhash::get_Bneighbor_substr2gap(string s)
{       
    int len = s.length();
    int k1 = klen - blen;

    readinfo* r = new readinfo(len, klen, k1);



    int next[4] = {-1, -1, -1, -1};
    int last[4] = {-1, -1, -1, -1};

    for(int i = 0; i < 4; i++)
    {
        edges[len-1][i] = -1;
        edgesv[0][i] = -1;
    }

    for(int i = len - 1; i > 0; i--)
    {
        int c = dict[s[i]];
        next[c] = i;

        for(int j = 0; j < 4; j++)
            edges[i-1][j] = next[j];
    }

    for(int i = 0; i < len - 1; i++)
    {
        int c = dict[s[i]];
        last[c] = i;

        for(int j = 0; j < 4; j++)
            edgesv[i+1][j] = last[j];
    }

    int num = len - klen + 1;
    int mnum = len - mlen + 1;

    if(num <= 0)
        return r;
    vector<int> indices;

    int x = 0;
    for(int i = 0; i < mlen - 1; i++)
        x = (x<<2) + dict[s[i]];

    for(int i = 0; i < mnum; i++)
    {
        x = x % (1<<(2*mlen-2));
        x = (x<<2) + dict[s[i + mlen - 1]];

        for(hashtable* j = heads[x]; j != NULL; j = j->next)
        {
            int id = j->rank;
            int pos = j->pos;
            int st = i, en = i + mlen - 1;

            int nbr = ranks[id]; 
            string c = bestb[id];

            int stext = pos - 1;
            int enext = pos + mlen;

            int exmatchst = i;
            int exmatchen = i + mlen - 1;

            indices.clear();
            for(int k = st; k <= en; k++)
                indices.push_back(k);

            for( ;stext >= 0; stext--)
            {
                int y = (nbr>>(2 * (k1 - stext - 1))) & 3;
                if(st != edgesv[st][y] + 1) 
                    break;

                exmatchst--;
                st = edgesv[st][y];
                indices.push_back(st);

                if(st == -1)
                    break;
            }

            for(; enext < k1; enext++)
            {
                int y = (nbr>>(2 * (k1 - enext - 1))) & 3;
                if(en != edges[en][y] - 1) 
                    break;

                en = edges[en][y];
                exmatchen++;
                indices.push_back(en);

                if(en == -1)
                    break;
            }

            if(st < 0 || en < 0)
                continue;


            if(stext >= 0 && enext < k1) //1 gap each side
            {
                st--;

                if(st >= 0)
                {
                    vector<int> pi = suffix_function(c.substr(0, stext + 1));
                    int max_len = max(0, st - (blen + stext));

                    int k = stext;
                    while(st >= max_len)
                    {
                        if(c[k] == s[st])
                        {
                            k--;
                            st--;

                            if(k == -1)
                            {
                                st++;
                                break;
                            }
                        }
                        else
                            if(k == stext)
                                st--;
                            else
                                k = stext - pi[stext - k - 1];
                    }

                    if(k != -1)
                        st = -1;                    
                    else
                        for(int kk = st; kk <= st + stext; kk++)
                            indices.push_back(kk);
                }

                en++;
                if(en < len)
                {
                    vector<int> pi = prefix_function(c.substr(enext));
                    int max_len = min(len - 1, en + (blen + k1 - enext) - 1);

                    int k = enext;
                    while(en <= max_len)
                    {
                        if(c[k] == s[en])
                        {
                            k++;
                            en++;

                            if(k == k1)
                            {
                                en--;
                                break;
                            }
                        }
                        else
                            if(k == enext)
                                en++;
                            else
                                k = enext + pi[k - enext - 1];
                    }

                    if(k != k1)
                        en = -1;
                    else
                        for(int kk = en  - (k1 - enext) + 1; kk <= en; kk++)
                            indices.push_back(kk);


                }
                else
                    en = -1;
            }

            else //2 gaps one side
            {
                if(stext >= 0)
                {               
                    st--;

                    if(st >= 0)
                    {
                        string ptr = c.substr(0, stext + 1);


                        int tst = max(0, st - (blen + stext));
                        int ten = st;

                        segment_gap mpos = prematch(s.substr(tst, ten - tst + 1), ptr);


                        if(mpos.st1 == -1)
                            st = -1;
                        else
                            st = tst + mpos.st1;

                        if(mpos.st1 != -1)
                            for(int k = mpos.st1; k <= mpos.en1; k++)
                                indices.push_back(tst + k);
                        if(mpos.st2 != -1)
                            for(int k = mpos.st2; k <= mpos.en2; k++)
                                indices.push_back(tst + k);
                    }
                }

                if(enext < k1 && en != -1)
                {
                    en++;
                    if(en < len)
                    {
                        string ptr = c.substr(enext);

                        int ten = min(len - 1, en + (blen + k1 - enext) - 1);
                        int tst = en;

                        segment_gap mpos = sufmatch(s.substr(tst, ten - tst + 1), ptr);

                        if(mpos.en1 == -1)
                            en = -1;
                        else
                        {   
                            if(mpos.st2 == -1)
                            {
                                en = tst + mpos.en1;
                            }
                            else
                            {
                                en = tst + mpos.en2;
                            }


                            if(mpos.st1 != -1)
                                for(int k = mpos.st1; k <= mpos.en1; k++)
                                    indices.push_back(tst + k);
                            if(mpos.st2 != -1)
                                for(int k = mpos.st2; k <= mpos.en2; k++)
                                    indices.push_back(tst + k);
                        }
                    }
                    else
                        en = -1;
                }
            }


            if(st < 0 || en < 0 || en - st + 1 > klen)
                continue;


            bnbr tmp;
            tmp.st = st;
            tmp.en = en;
            tmp.est = i;
            tmp.een = i + mlen - 1;

            for(int index:indices)
                tmp.indices.push_back(index);
            
            tmp.id = id;

            if(r->slist[st][en-st+1-k1].id == -1 || r->slist[st][en-st+1-k1].id > id)
                r->slist[st][en-st+1-k1] = tmp;               
                
        }
    }

    r->getminbnlist();
    //concat(r, s);
    return r;
}


// bool compare_match(matchpoint a, matchpoint b)
// {
//     return a.st1 < b.st1 || (a.st1 == b.st1 && a.st2 > b.st2);
// }



// vector<matchpoint> kbhash::chaining(vector<matchpoint> matches)
// {
//     vector<matchpoint> ret;

//     int m = matches.size();
//     for(int i = 0; i < m; i++)
//     {
//         g[i] = -1;
//         dp[i] = 1;
//     }   

//     int ans = 0, pos;

//     for(int i = 0; i < m; i++)
//     {
//         for(int j = 0; j < i; j++)
//             if(matches[i].st2 > matches[j].st2
//                 && dp[j] + 1  > dp[i])
//             {
//                 dp[i] = dp[j] + 1;
//                 g[i] = j;
//             }

//         if(dp[i] > ans)
//         {
//             ans = dp[i];
//             pos = i;
//         }
//     }

//     //cout<<pos<<" "<<ans<<endl;
//     int en1 = matches[pos].en1;
//     int en2 = matches[pos].en2;
//     int noc = 0;
//     double avgdis = 0;


//     while(g[pos] != -1)
//     {

//         ret.push_back(matches[pos]);
//         // for(int i = matches[pos].st1; i <= matches[pos].en1; i++)
//         //     mcv1[i] = '0' + (noc % 10);

//         // for(int i = matches[pos].st2; i <= matches[pos].en2; i++)
//         //     mcv2[i] = '0' + (noc % 10);

//         // for(int i : matches[pos].b1.indices)
//         //     scv1[i] = '0' + (noc % 10);
//         // for(int i : matches[pos].b2.indices)
//         //     scv2[i] = '0' + (noc % 10);
//         avgdis += matches[pos].st1 - matches[pos].st2;
//         cout<<matches[pos].st1<<" "<<matches[pos].st2<<" "<<matches[pos].st1 - matches[pos].st2<<endl;
//         noc++;
//         //cout<<matches[pos].st1<<" "<<matches[pos].st2<<endl;
//         pos = g[pos];
//     }
//     avgdis /= ret.size();
//     //printf("%.2lf\n", avgdis);
//     // for(int i = matches[pos].st1; i <= matches[pos].en1; i++)
//     //     mcv1[i] = '0' + (noc % 10);

//     // for(int i = matches[pos].st2; i <= matches[pos].en2; i++)
//     //     mcv2[i] = '0' + (noc % 10);

//     //     for(int i : matches[pos].b1.indices)
//     //         scv1[i] = '0' + (noc % 10);
//     //     for(int i : matches[pos].b2.indices)
//     //         scv2[i] = '0' + (noc % 10);
//     noc++;

//     int st1 = matches[pos].st1;
//     int st2 = matches[pos].st2;

//     //printf("(%d, %d), %d\n", st1, en1, len1);
//     //printf("(%d, %d), %d\n", st2, en2, len2);

//     // int len1 = 0, len2 = 0;
//     // for(int i = 0; i <= en1; i++)
//     // {
//     //     if(scv1[i] != '-')
//     //         len1++;
//     //     if(mcv1[i] != '-')
//     //         totallen++;
//     // }

//     // for(int i = 0; i <= en2; i++)
//     // {
//     //     if(scv2[i] != '-')
//     //         len2++;
//     //     if(mcv2[i] != '-')
//     //         totallen++;
//     // }



//     //totalnum += len1 + len2;
//     //totalnoc += noc;
//     //cout<<len1 + len2<<endl;
//     //cout<<ans<<endl;
//     int nst = 0, nen = ret.size() - 1;
//     while(ret[nst].st1 - ret[nst+1].st1 > 20 || ret[nst].st2 - ret[nst+1].st2 > 20)
//     {
//         // cout<<ret[nst].st1<<" "<<ret[nst].st2<<endl;
//         // cout<<ret[nst+1].st1<<" "<<ret[nst+1].st2<<endl;
//         // cout<<endl;
//         nst++;
//     }

//     while(ret[nen-1].st1 - ret[nen].st1 > 20 || ret[nen-1].st2 - ret[nen].st2 > 20)
//     {
//         // cout<<ret[nen -1].st1<<" "<<ret[nen -1].st2<<endl;
//         // cout<<ret[nen].st1<<" "<<ret[nen].st2<<endl;
//         // cout<<endl;
//         nen--;
//     }

//     vector<matchpoint> nret;
//     for(int i = nen; i >= nst; i--)
//         nret.push_back(ret[i]);

//     //cout<<ret.size()<<" "<<nret.size()<<endl;
//     return nret;
// }

// vector<matchpoint> kbhash::match(readinfo* r1, readinfo* r2)
// {
//     int size = r1->len - klen + 1;
    

//     set<int> s1, s2;
//     set<int> ms1, ms2;
//     set<int> ss1, ss2;

//     for(auto b: r1->minbnlist)
//     {
//         if(s1.find(b.id) != s1.end())
//             ms1.insert(b.id);
//         else
//             s1.insert(b.id);
//     }
//     set_difference(s1.begin(), s1.end(), ms1.begin(), ms1.end(), inserter(ss1 , ss1.begin()));

//     for(auto b: r2->minbnlist)
//     {
//         if(s2.find(b.id) != s2.end())
//             ms2.insert(b.id);
//         else
//             s2.insert(b.id);
//     }
//     set_difference(s2.begin(), s2.end(), ms2.begin(), ms2.end(), inserter(ss2 , ss2.begin()));


//     set<int> is;
//     set_intersection(ss1.begin(), ss1.end(), ss2.begin(), ss2.end(), inserter(is , is.begin()));

//     //cout<<s1.size()<<" "<<s2.size()<<" "<<is.size()<<endl;
//     totalset += s1.size() + s2.size();
//     totalinter += is.size();
    
//     //get all the matches
//     vector<matchpoint> matches;
//     map<int,int> bmap;

//     for(int b: is)
//     {
//         matchpoint tmp;
//         tmp.id = b;
//         bmap.insert(make_pair(b, matches.size()));
//         matches.push_back(tmp);
//     }

//     for(auto b: r1->minbnlist)
//         if(bmap.find(b.id) != bmap.end())
//         {
//             int pos = bmap[b.id];
//             matches[pos].st1 = b.vst;
//             matches[pos].en1 = b.ven;
//             matches[pos].b1 = b;
//         }

//     for(auto b: r2->minbnlist)
//         if(bmap.find(b.id) != bmap.end())
//         {
//             int pos = bmap[b.id];
//             matches[pos].st2 = b.vst;
//             matches[pos].en2 = b.ven;
//             matches[pos].b2 = b;
//         }

//     //sort all the matches then merge
    
//     // sort(matches.begin(), matches.end(), compare_point);

//     // //merge overlapped matches
//     // matchpoint tmp;
//     // tmp.id = -1;
//     // vector<matchpoint> merged;

//     // for(auto b: matches)
//     // {
//     //  //cout<<b.b1.st<<" "<<b.b1.en<<" "<<b.b2.st<<" "<<b.b2.en<<" "<<b.id<<endl;
//     //  if(tmp.id == -1)
//     //  {
//     //      tmp = b;
//     //      tmp.st1 = b.b1.st;
//     //      tmp.en1 = b.b1.en;
//     //      tmp.st2 = b.b2.st;
//     //      tmp.en2 = b.b2.en;
//     //      continue;
//     //  }   

//     //  if(b.b1.st <= tmp.en1 && b.b2.st > tmp.st2 && b.b2.st <= tmp.en2)
//     //  {
//     //      tmp.en1 = b.b1.en;
//     //      tmp.en2 = b.b2.en;
//     //      tmp.num++;
//     //  }
//     //  else
//     //  {
//     //      merged.push_back(tmp);
//     //      tmp = b;
//     //      tmp.st1 = b.b1.st;
//     //      tmp.en1 = b.b1.en;
//     //      tmp.st2 = b.b2.st;
//     //      tmp.en2 = b.b2.en;
//     //  }
//     // }

//     // if(tmp.id != -1)
//     //  merged.push_back(tmp);

//     //Multiple Matches splited into single matches
//     bmap.clear();
//     vector<multiple> multimatch;
//     for(int b: ms1)
//         if(s2.find(b) != s2.end())
//         {
//             multiple tmp;
//             tmp.id = b;
//             bmap.insert(make_pair(b, multimatch.size()));
//             multimatch.push_back(tmp);
//         }
//     for(int b: ms2)
//         if(s1.find(b) != s1.end() && bmap.find(b) == bmap.end())
//         {
//             multiple tmp;
//             tmp.id = b;
//             bmap.insert(make_pair(b, multimatch.size()));
//             multimatch.push_back(tmp);
//         }

//     for(auto b: r1->minbnlist)
//         if(bmap.find(b.id) != bmap.end())
//         {
//             int pos = bmap[b.id];
//             multimatch[pos].b1.push_back(b);
//         }

//     for(auto b: r2->minbnlist)
//         if(bmap.find(b.id) != bmap.end())
//         {
//             int pos = bmap[b.id];
//             multimatch[pos].b2.push_back(b);
//         }

//     for(auto x: multimatch)
//         for(auto b1: x.b1)
//             for(auto b2: x.b2)
//             {
//                 matchpoint tmp;
//                 tmp.id = x.id;
//                 tmp.st1 = b1.vst;
//                 tmp.en1 = b1.ven;
//                 tmp.st2 = b2.vst;
//                 tmp.en2 = b2.ven;
//                 tmp.b1 = b1;
//                 tmp.b2 = b2;
//                 matches.push_back(tmp);
//             }

//     totalanc += matches.size();
//     sort(matches.begin(), matches.end(), compare_match);
//     //cout<<matches.size()<<endl;
//     //if(matches.size() > 0)
//     return chaining(matches);
// }