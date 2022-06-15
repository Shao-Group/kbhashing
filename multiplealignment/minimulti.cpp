#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <climits>
#include <map>
#include "dpchain.h"


using namespace std;

map<char,int> dict;

uint64_t get_hash(string s, int pos, int len)
{
	uint64_t ret = 0;
	for(int i = 0; i < len; i++)
		ret = ret * 4 + dict[s[pos + i]];

    uint64_t m = 0xFFFFFFFFFFFFFFFF;
    //cout<<m<<endl;

    ret = (~ret + (ret << 21)) & m;
    ret = ret ^ ret >> 24;
    ret = (ret + (ret<<3) + (ret<<8)) & m;
    ret = ret ^ ret >> 14;
    ret = (ret + (ret<<2) + (ret<<4)) & m;
    ret = ret ^ ret >> 28;
    ret = (ret + (ret<<31)) & m;

	return ret;
}

vector<seed> get_minimizers(string s, int k, int w)
{
	int len = s.length();
	vector<seed> ret;
	vector<uint64_t> mins;

	for(int i = 0; i <= len - k; i++)
	{
		uint64_t v = get_hash(s, i, k);
		mins.push_back(v);
		//cout<<i<<" "<<v<<" "<<s.substr(i, k)<<endl;
	}

	for(int i = 0; i <= len - k - w + 1; i++)
	{
		uint64_t minimum = 0xFFFFFFFFFFFFFFFF;
		int pos = i;

		for(int j = i; j < i + w; j++)
			if(mins[j] < minimum)
			{
				minimum = mins[j];
				pos = j;
			}

		if(ret.empty() || pos != ret.back().st || minimum != ret.back().hashval)
		{
			seed tmp;
			tmp.st = pos;
			tmp.en = pos + k - 1;
			tmp.hashval = minimum;

			//cout<<pos<<" "<<minimum<<endl;
			for(int j = pos; j < pos + k; j++)
				tmp.index.push_back(j);

			// cout<<pos<<" "<<tmp.hashval<<endl;
			ret.push_back(tmp);
		}

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



void minimizer_match(int i1, int i2, int k, int w, vector<int> align)
{
	int trueset = 0;

	for(int i: align)
		trueset += (i != -1);

	clock_t t0, t1, t2;
    t0 = clock();
	vector<matchpoint> matches = chain.getmatches(seeds[i1], seeds[i2], ancs[i1], ancs[i2], multiancs[i1], multiancs[i2], singleancs[i1], singleancs[i2]);

	totalanc += matches.size();
	
	result ret;
	// cout<<matches.size()<<endl;
	// for(auto a:matches)
	// 	cout<<a.hashval<<" "<<a.ac1.st<<" "<<a.ac1.en<<" "<<a.ac2.st<<" "<<a.ac2.en<<endl;

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
	int w = stoi(argv[2]);

	string s, t;
	int n = 0;

    dict['A'] = 0;
    dict['C'] = 1;
    dict['G'] = 2;
    dict['T'] = 3;

	int tmp, len;    
	clock_t t0, t1, t2;
    t0 = clock();

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
		seeds.push_back(get_minimizers(t, k, w));
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

			minimizer_match(j, i, k, w, align);
		}
	}

    t2 = clock();
	printf("%.2lf\n", double(t2 - t1)/CLOCKS_PER_SEC);
	int m = n * (n-1) / 2;
	printf("%d, %d, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf\n", k, w, totalanc / m, totalacc / m, tp1 / m, rc1 / m, bp1 / m, tp2 / m, rc2 / m, bp2 / m, totalnum / m, 1 - totalnum / m, totalil / m, abp / m);
	//printf("%.2lf, %.2lf, %.2lf, %.2lf, %.2lf\n", wm / m, manc / m, avgac / m, segnum / m, avgseg / m);
	printf("%.2lf, %.2lf, %.2lf\n", time0 / CLOCKS_PER_SEC, time1  / CLOCKS_PER_SEC, time2  / CLOCKS_PER_SEC);

	return 0;
}