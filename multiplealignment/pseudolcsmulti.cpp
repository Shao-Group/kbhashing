#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <climits>
#include <map>
#include "dpchain.h"

random_device rd;

uint64_t f[50][50];
int g[50][50];
map<char, uint64_t> dict[50];
vector<int> pos;

void init(int len, int blen)
{
	pos.clear();
	for(int i = 0; i < blen; i++)
		pos.push_back(i);

	unsigned seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
	shuffle(pos.begin(), pos.end(), default_random_engine(seed));

	vector<char> ab;
	ab.push_back('A');
	ab.push_back('C');
	ab.push_back('G');
	ab.push_back('T');


	for(int i = 0; i < blen; i++)
	{
		seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
		shuffle(ab.begin(), ab.end(), default_random_engine(seed));
		for(int j = 0; j < 4; j++)
			dict[i][ab[j]] = j;
	}
}

matchpoint get_lac(matchpoint a)
{
	int ansst = 1, ans = 1;
	int st = 0, tl = 1;
	int len = a.ac1->index.size();
	for(int i = 1; i < len; i++)
	{
		if(a.ac1->index[i] == a.ac1->index[i-1] - 1 && a.ac2->index[i] == a.ac2->index[i-1] - 1)
			tl++;
		else
		{
			st = i;
			tl = 1;
		}

		if(tl > ans)
		{
			ans = tl;
			ansst = st;
		}
	}

	matchpoint tmp;
	tmp.hashval = a.hashval;
	seed* tmp1 = new seed(), *tmp2 = new seed();
	for(int i = 0; i < ans; i++)
	{
		tmp1->index.push_back(a.ac1->index[ansst+i]);
		tmp2->index.push_back(a.ac2->index[ansst+i]);
	}
	tmp1->st = tmp1->index.back();
	tmp1->en = tmp1->index.front();
	tmp2->st = tmp2->index.back();
	tmp2->en = tmp2->index.front();

	tmp.ac1 = tmp1;
	tmp.ac2 = tmp2;

	return tmp;
}

vector<seed> get_bmers(string s, int k, int b)
{
	int del = k - b;
	vector<seed> ret;

	f[0][0] = 0;
	for(int i = 1; i <= k; i++)
	{
		f[i][0] = 0;
		for(int j = 1; j <= b; j++)
			f[i][j] = LLONG_MAX>>1;
	}

	int len = s.length();

	for(int st = 0; st + k <= len; st++)
	{
		f[1][1] = dict[0][s[st]] << (2*pos[0]);
		g[1][1] = 0;

		for(int i = 2; i <= k; i++)
		{
			int minj = max(1, i - del);
			int maxj = min(i, b);

			for(int j = minj; j <= maxj; j++)
			{
				f[i][j] = f[i-1][j];
				g[i][j] = g[i-1][j];

				if(f[i-1][j-1] + (dict[j-1][s[st + i - 1]] << (2*pos[j-1])) < f[i][j])
				{
					f[i][j] = f[i-1][j-1] + (dict[j-1][s[st + i - 1]] << (2*pos[j-1]));
					g[i][j] = i-1;
				}
			}
		}

		if(ret.size() > 0 && f[k][b] == ret.back().hashval && ret.back().index[0] == g[k][b] + st)
			continue;

		seed tmp;
		tmp.hashval = f[k][b];
		int x = k;
		int y = b;


		while(y > 0)
		{
			tmp.index.push_back(g[x][y] + st);
			x = g[x][y];
			y--;
		}

		tmp.st = tmp.index.back();
		tmp.en = tmp.index.front();

		// for(auto c: tmp.index)
		// 	cout<<c<<" ";
		// cout<<endl;
		//tmp = get_lac(tmp);
		// for(auto c: tmp.index)
		// 	cout<<c<<" ";
		// cout<<endl;
		// cout<<endl;
		ret.push_back(tmp);
	}

	return ret;
}

double time0 = 0, time1 = 0, time2 = 0;
double totalanc = 0, manc = 0, totalacc = 0, totalnum = 0, totalil = 0, abp = 0, wm = 0;
double tp1 = 0, tp2 = 0 , rc1 = 0, rc2 = 0;
double bp1 = 0, bp2 = 0;
double avgseg = 0, avgac = 0, segnum = 0;
dpchain chain;


vector<string> reads;
vector<vector<seed>> seeds;

vector<set<uint64_t>> ancs;
vector<set<uint64_t>> multiancs;
vector<set<uint64_t>> singleancs;


void pseudo_match(int i1, int i2, int k, int blen, vector<int> align)
{
	int trueset = 0;

	for(int i: align)
		trueset += (i != -1);

	clock_t t0, t1, t2;
    t0 = clock();
	vector<matchpoint> matches = chain.getmatches(seeds[i1], seeds[i2], ancs[i1], ancs[i2], multiancs[i1], multiancs[i2], singleancs[i1], singleancs[i2]);
	// for(auto a:matches)
	// 	cout<<a.hashval<<" "<<a.ac1->st<<" "<<a.ac1->en<<" "<<a.ac2->st<<" "<<a.ac2->en<<endl;

	for(int i = 0; i < matches.size(); i++)
		matches[i] = get_lac(matches[i]);

	totalanc += matches.size();
	
	result ret;
	// cout<<matches.size()<<endl;

    t1 = clock();

	if(matches.size() > 0)
		ret = chain.alignment(matches, align, reads[i1].length(), reads[i2].length(), k, reads[i1], reads[i2]);
	
	
	totalacc += (ret.tp1 + ret.tp2)/trueset;

	if(trueset != ret.ts)
		tp1 += ret.tp1/(trueset - ret.ts);
	if(ret.ts != 0)
		tp2 += ret.tp2/ret.ts;

	if(ret.pd1 != 0)
		rc1 += ret.tp1/ret.pd1;
	if(ret.pd2 != 0)
		rc2 += ret.tp2/ret.pd2;

	bp1 += ret.tp1 / reads[i1].length();
	bp2 += ret.tp2 / reads[i1].length();

	manc += ret.num;
	wm += ret.wm;

	totalil += ret.island;
	totalnum += ret.rt;
	abp += ret.abproduct;

	avgseg += ret.avgseg;
	avgac += ret.avgac;
	segnum += ret.segnum;

	time0 += double(t1 - t0);
	time1 += ret.t1;
	time2 += ret.t2;
}

vector<vector<int>> aligns;
int refalign[12000];
const int reflen = 10000;

int main(int argc, const char * argv[])
{
	int k = stoi(argv[1]);
	int b = stoi(argv[2]);

	string s, t;
	int n = 0;

	int tmp, len;    
	clock_t t0, t1, t2;
    t0 = clock();

	init(k, b);

	cin>>s;

	while(cin>>t)
	{
    	vector<int> align;

		len = t.length();
		for(int i = 0; i < len; i++)
		{
			cin>>tmp;
			align.push_back(tmp);
		}

		reads.push_back(t);
		seeds.push_back(get_bmers(t, k, b));
		aligns.push_back(align);

		set<uint64_t> t1, ms1, ss1;

		for(auto b: seeds[n])
		{
			if(t1.find(b.hashval) != t1.end())
				ms1.insert(b.hashval);
			else
				t1.insert(b.hashval);
		}
		set_difference(t1.begin(), t1.end(), ms1.begin(), ms1.end(), inserter(ss1, ss1.begin()));

		ancs.push_back(t1);
		multiancs.push_back(ms1);
		singleancs.push_back(ss1);

		n++;
	}


    t1 = clock();
    printf("%.2lf\n", double(t1 - t0)/CLOCKS_PER_SEC);

	for(int i = 0; i < n; i++)
	{

		len = reads[i].length();

		for(int l = 0; l < reflen; l++)
			refalign[l] = -1;
		for(int l = 0; l < len; l++)
			if(aligns[i][l] != -1)
				refalign[aligns[i][l]] = l;

		for(int j = i+1; j < n; j++)
		{
			len = reads[j].length();
			vector<int> align;

			for(int l = 0; l < len; l++)
				if(aligns[j][l] != -1)
					align.push_back(refalign[aligns[j][l]]);
				else
					align.push_back(-1);

			pseudo_match(j, i, k, b, align);
		}
	}

    t2 = clock();
	printf("%.2lf\n", double(t2 - t1)/CLOCKS_PER_SEC);
	int m = n * (n-1) / 2;
	printf("%d, %d, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf\n", k, b, totalanc / m, totalacc / m, tp1 / m, rc1 / m, bp1 / m, tp2 / m, rc2 / m, bp2 / m, totalnum / m, 1 - totalnum / m, totalil / m, abp / m);
	//printf("%.2lf, %.2lf, %.2lf, %.2lf, %.2lf\n", wm / m, manc / m, avgac / m, segnum / m, avgseg / m);
	printf("%.2lf, %.2lf, %.2lf\n", time0 / CLOCKS_PER_SEC, time1  / CLOCKS_PER_SEC, time2 / CLOCKS_PER_SEC);

	return 0;
}