#include "dpchain.h"

using namespace std;


int dpchain::binary_search(int s, int t, int v)
{
	while(s < t)
	{
		int mid = (s+t)>>1;
		if(dp[mid] < v)
			s = mid + 1;
		else 
			t = mid;
	}

	return s;
}

pair<double, double> dpchain::chaining(vector<matchpoint> matches, vector<int> align, int lens, int lent, int blen)
{	
    int m = matches.size();
	for(int i = 0; i < m; i++)
		g[i] = -1;

	int ans = 1, pos;
	dp[1] = matches[0].ac2->st;
	idx[1] = 0;

	for(int i = 1; i < m; i++)
	{
		if(matches[i].ac2->st > dp[ans])
		{
			g[i] = idx[ans];
			ans++;
			dp[ans] = matches[i].ac2->st;
			idx[ans] = i;
		}

		else
		{
			int j = binary_search(1, ans, matches[i].ac2->st);
		
			if(j > 1)
				g[i] = idx[j-1];

			if(matches[i].ac2->st < dp[j])
			{
				dp[j] = matches[i].ac2->st;
				idx[j] = i;
			}
		}
	}

	pos = idx[ans];

 	vector<int> salign(lens, 0);
 	vector<int> talign(lens, 0);


    while(pos != -1)
    {
        for(int i:matches[pos].ac1->index)
            scv1[i] = '1';

        for(int i:matches[pos].ac2->index)
            scv2[i] = '1';

        for(int i = 0; i < blen; i++)
        {
        	int p1 = matches[pos].ac1->index[i];
        	int p2 = matches[pos].ac2->index[i];

        	//cout<<p1<<" "<<p2<<endl;
        	if(align[p1] == p2)
        		talign[p1]++;
        	salign[p1]++;
        }


        pos = g[pos];
    }

    int len1 = 0, len2 = 0;

    for(int i = 0; i < lens; i++)
        if(scv1[i] != '-')
            len1++;
    for(int i = 0; i < lent; i++)
        if(scv2[i] != '-')
            len2++;

    double tpalign = 0;
    for(int i = 0; i < lens; i++)
    	if(align[i] != -1 && salign[i] > 0)
    		tpalign += talign[i]/salign[i];

    return make_pair(tpalign/double(lens), double(len1)/(lens) + double(len2)/(lent));
}

vector<pair<int,int>> dpchain::GlobalAlignment(string kmer1, string kmer2, int st1, int en1, int st2, int en2)
{
    int ans = 0, posx, posy;

    int len1 = en1 - st1 + 1;
    int len2 = en2 - st2 + 1;

    f[0][0] = 0;

    for(int i = 1; i <= len1; i++)
    	f[i][0] = i;

    for(int i = 1; i <= len2; i++)
    	f[0][i] = i;

    for(int i = 1; i <= len1; i++)
        for(int j = 1; j <= len2; j++)
        {
            if(kmer1[st1+i-1] == kmer2[st2+j-1])
                f[i][j] = f[i-1][j-1] + 1;
            else
                f[i][j] = f[i-1][j-1] - 1;
        
            f[i][j] = max(f[i][j], max(f[i-1][j] - 1, f[i][j-1] - 1));

        }
    ans = f[len1][len2];
    int x = len1, y = len2;

    vector<pair<int,int>> v;

    while(x > 0 &&  y > 0)
    {
        if(kmer1[st1+x-1] == kmer2[st2+y-1] && f[x][y] == f[x-1][y-1] + 1)
        {
        	v.push_back(make_pair(st1 + x - 1,st2 + y - 1));
            x--;
            y--;
        }
        else if(f[x][y] == f[x-1][y-1] - 1)
        {
            x--;
            y--;
        }
        else if(f[x][y] == f[x-1][y] - 1)
        {
            x--;
        }
        else if(f[x][y] == f[x][y-1] - 1)
        {
            y--;
        }

    }

    return v;
}

result dpchain::alignment(vector<matchpoint> &matches, vector<int> &align, int lens, int lent, int blen, string s, string t)
{	
	clock_t t0, t1, t2, t3;
    t0 = clock();
    int m = matches.size();
	for(int i = 0; i < m; i++)
		g[i] = -1;

	int ans = 1, pos;
	dp[1] = matches[0].ac2->st;
	idx[1] = 0;

	for(int i = 1; i < m; i++)
	{
		if(matches[i].ac2->st > dp[ans])
		{
			g[i] = idx[ans];
			ans++;
			dp[ans] = matches[i].ac2->st;
			idx[ans] = i;
		}

		else
		{
			int j = binary_search(1, ans, matches[i].ac2->st);
		
			if(j > 1)
				g[i] = idx[j-1];

			if(matches[i].ac2->st < dp[j])
			{
				dp[j] = matches[i].ac2->st;
				idx[j] = i;
			}
		}
	}
	
    t1 = clock();

	pos = idx[ans];

 	vector<interval> intervals;
 	vector<int> salign(lens, 0);
 	vector<int> talign(lens, 0);
 	int wm = 0, num = 0;

 	if(matches[pos].ac1->en < lens - 1 && matches[pos].ac2->en < lent - 1)
 	{
 		interval inv;
 		inv.st1 = matches[pos].ac1->en + 1;
 		inv.en1 = lens - 1;
 		inv.st2 = matches[pos].ac2->en + 1;
 		inv.en2 = lent - 1;
 		intervals.push_back(inv);
 	}

 	double avgac = 0;
 	double avgseg = 0;
 	int segnum = 0;

    while(pos != -1)
    {
    	//cout<<matches[pos].ac1->st<<" "<<matches[pos].ac1->en<<" "<<matches[pos].ac2->st<<" "<<matches[pos].ac2->en<<endl;
    	// cout<<s.substr(matches[pos].ac1->st, matches[pos].ac1->en - matches[pos].ac1->st + 1)<<endl;
    	// cout<<t.substr(matches[pos].ac2->st, matches[pos].ac2->en - matches[pos].ac2->st + 1)<<endl;
    	avgac += matches[pos].ac1->en - matches[pos].ac1->st + 1;

    	int tmpnum = 0;
    	num++;

        for(int i = 0; i < matches[pos].ac1->index.size(); i++)
        {
        	int p1 = matches[pos].ac1->index[i];
        	int p2 = matches[pos].ac2->index[i];

        	if(align[p1] == p2)
        	{
        		talign[p1]++;
        	}

        	if(abs(align[p1] - p2) <= 10)
        		tmpnum++;

        	salign[p1]++;
        }


        if(tmpnum == 0)
        {
    		// cout<<matches[pos].ac1->st<<" "<<matches[pos].ac1->en<<" "<<matches[pos].ac2->st<<" "<<matches[pos].ac2->en<<endl;
	     //    for(int i = 0; i < matches[pos].ac1->index.size(); i++)
	     //    {
      //   	int p1 = matches[pos].ac1->index[i];
      //   	int p2 = matches[pos].ac2->index[i];
      //   	cout<<p1<<" "<<align[p1]<<" "<<p2<<endl;
      //   }
        	wm += 1;
        }

        if(g[pos] != -1)
        {
        	if(matches[g[pos]].ac1->en < matches[pos].ac1->st - 1 && matches[g[pos]].ac2->en < matches[pos].ac2->st - 1)
		 	{
		 		interval inv;
		 		inv.st1 = matches[g[pos]].ac1->en + 1;
		 		inv.en1 = matches[pos].ac1->st - 1;
		 		inv.st2 = matches[g[pos]].ac2->en + 1;
		 		inv.en2 = matches[pos].ac2->st - 1;
		 		intervals.push_back(inv);
		 	}
        }
        else
        {
        	if(matches[pos].ac1->st > 0 && matches[pos].ac2->st > 0)
        	{
		 		interval inv;
		 		inv.st1 = 0;
		 		inv.en1 = matches[pos].ac1->st - 1;
		 		inv.st2 = 0;
		 		inv.en2 = matches[pos].ac2->st - 1;
		 		intervals.push_back(inv);

        	}
        }

        pos = g[pos];
    }


    double tpalign1 = 0, tpalign2 = 0;
    double island = 0;
    int islandnum = 0;
    double interval_len = 0;
    double abproduct = 0;
    int trueset2 = 0;
    int preditcion1 = 0, preditcion2 = 0;

    int last = lens - 1;

    for(auto i: intervals)
    {
    	vector<pair<int,int>> ap = GlobalAlignment(s, t, i.st1, i.en1, i.st2, i.en2);

    	//cout<<i.st1<<" "<<i.en1<<" "<<i.st2<<" "<<i.en2<<endl;
    	// cout<<island<<endl;
    	island += (i.en1 - i.st1 + 1) * (i.en1 - i.st1 + 1) + (i.en2 - i.st2 + 1) * (i.en2 - i.st2 + 1);
    	interval_len += (i.en1 - i.st1 + 1) + (i.en2 - i.st2 + 1);
    	abproduct += (i.en1 - i.st1 + 1) * (i.en2 - i.st2 + 1);
    	islandnum += 2;

    	preditcion1 += ap.size();

    	for(auto pairs:ap)
    	{
    		if(align[pairs.first] == pairs.second)
    			tpalign1 += 1;
    	}

    	avgseg += last - i.en1;
    	if(i.en1 < last)
    	{

    		segnum++;
    	}

    	for(int j = last; j > i.en1; j--)
	    {
	  		if(align[j] != -1)
	    		trueset2++;
	    	if(salign[j] > 0)
	    		preditcion2++;
	    	if(align[j] != -1 && salign[j] > 0 && talign[j] * 2 > salign[j])
	    		tpalign2 += 1;
	    }

	    last = i.st1 - 1;
    }

    if(last >= 0)
    {
    	segnum ++;
    	avgseg += last;
    	for(int j = last; j >= 0; j--)
	    {
	  		if(align[j] != -1)
	    		trueset2++;
	    	if(salign[j] > 0)
	    		preditcion2++;
	    	if(align[j] != -1 && salign[j] > 0 && talign[j] * 2 > salign[j])
	    		tpalign2 += 1;
	    }
	}

    t2 = clock();
	//printf("%.2lf\t%.2lf\t%.2lf\n", double(t1 - t0)/CLOCKS_PER_SEC, double(t2 - t1)/CLOCKS_PER_SEC, double(t3 - t2)/CLOCKS_PER_SEC);
    result ret;
    
    ret.tp1 = tpalign1;
    ret.tp2 = tpalign2;
    ret.ts = trueset2;
    ret.pd1 = preditcion1;
    ret.pd2 = preditcion2;
    ret.rt = interval_len/(lens + lent);
    ret.island = island/islandnum;
    ret.abproduct = abproduct;
    ret.wm = wm;
    ret.num = num;
    ret.avgac = avgac/num;
    ret.avgseg = avgseg/segnum;
    ret.segnum = segnum;

    ret.t1 = double(t1 - t0);
    ret.t2 = double(t2 - t1);
    return ret;
}

vector<matchpoint> dpchain::getmatches(vector<seed> &s1, vector<seed> &s2, set<uint64_t> &t1, set<uint64_t> &t2, set<uint64_t> &ms1, set<uint64_t> &ms2, set<uint64_t> &ss1, set<uint64_t> &ss2)
{
	set<uint64_t> is;
	set_intersection(ss1.begin(), ss1.end(),ss2.begin(), ss2.end(), inserter(is, is.begin()));

	vector<matchpoint> matches;
	unordered_map<uint64_t, int> bmap;

	for(uint64_t b: is)
	{
		matchpoint tmp;
		tmp.hashval = b;
		bmap.insert(make_pair(b, matches.size()));
		matches.push_back(tmp);
	}

	for(int i = 0; i < s1.size(); i++)
		if(bmap.find(s1[i].hashval) != bmap.end())
		{
			int pos = bmap[s1[i].hashval];
			matches[pos].ac1 = &s1[i];
		}

	for(int i = 0; i < s2.size(); i++)
		if(bmap.find(s2[i].hashval) != bmap.end())
		{
			int pos = bmap[s2[i].hashval];
			matches[pos].ac2 = &s2[i];
		}

	//Multiple Matches splited into single matches
	bmap.clear();
	vector<multiple> multimatch;
	for(uint64_t b: ms1)
		if(t2.find(b) != t2.end())
		{
			multiple tmp;
			tmp.hashval = b;
			bmap.insert(make_pair(b, multimatch.size()));
			multimatch.push_back(tmp);
		}
	for(uint64_t b: ms2)
		if(t1.find(b) != t1.end() && bmap.find(b) == bmap.end())
		{
			multiple tmp;
			tmp.hashval = b;
			bmap.insert(make_pair(b, multimatch.size()));
			multimatch.push_back(tmp);
		}

	for(int i = 0; i < s1.size(); i++)
		if(bmap.find(s1[i].hashval) != bmap.end())
		{
			int pos = bmap[s1[i].hashval];
			multimatch[pos].ac1.push_back(&s1[i]);
		}

	for(int i = 0; i < s2.size(); i++)
		if(bmap.find(s2[i].hashval) != bmap.end())
		{
			int pos = bmap[s2[i].hashval];
			multimatch[pos].ac2.push_back(&s2[i]);
		}

	for(auto x: multimatch)
		for(int i = 0; i < x.ac1.size(); i++)
			for(int j = 0; j < x.ac2.size(); j++)
			{
				matchpoint tmp;
				tmp.ac1 = x.ac1[i];
				tmp.ac2 = x.ac2[j];
				tmp.hashval = x.hashval;

				matches.push_back(tmp);
			}

	sort(matches.begin(), matches.end(), compare_match);

	return matches;
}

vector<matchpoint> dpchain::getmatches2(vector<seed> s1, vector<seed>s2, set<uint64_t> t1, set<uint64_t> t2)
{
	set<uint64_t> is;
	set_intersection(t1.begin(), t1.end(),t2.begin(), t2.end(), inserter(is, is.begin()));

	vector<matchpoint> matches;
	vector<multiple> multimatch;
	map<uint64_t, int> bmap;

	//Multiple Matches splited into single matches
	
	for(uint64_t b: is)
	{
		multiple tmp;
		tmp.hashval = b;
		bmap.insert(make_pair(b, multimatch.size()));
		multimatch.push_back(tmp);
	}

	for(int i = 0; i < s1.size(); i++)
		if(bmap.find(s1[i].hashval) != bmap.end())
		{
			int pos = bmap[s1[i].hashval];
			multimatch[pos].ac1.push_back(&s1[i]);
		}

	for(int i = 0; i < s2.size(); i++)
		if(bmap.find(s2[i].hashval) != bmap.end())
		{
			int pos = bmap[s2[i].hashval];
			multimatch[pos].ac2.push_back(&s2[i]);
		}


	for(auto x: multimatch)
		for(int i = 0; i < x.ac1.size(); i++)
			for(int j = 0; j < x.ac2.size(); j++)
			{
				matchpoint tmp;
				tmp.ac1 = x.ac1[i];
				tmp.ac2 = x.ac2[j];
				tmp.hashval = x.hashval;

				matches.push_back(tmp);
			}

	sort(matches.begin(), matches.end(), compare_match);

	return matches;
}
