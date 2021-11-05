#include <iostream>
#include <random>
#include <map>
#include <fstream>

using namespace std;
random_device rd;

char nt[4] = {'A', 'C', 'G', 'T'};
map<char,int> dict;

string simulation(int len)
{
	string s = "";
	for(int i = 0; i < len; i++)
	{
		int x = rd() % 4;
		s = s + nt[x];
	}

	return s;
}

string c_mutation(string s, int step) //Fixed position mutation
{
	int i, len = s.length();
	string t;
	for(i = step - 1; i < len; i += step)
	{
		int x = rd() % 3;
		t = t + s.substr(i - step + 1, step - 1);

		if(x == 0)
		{
			char c = nt[(dict[s[i]] + (rd() % 3) + 1) % 4];
			t = t + c;
		}
		else if(x == 1)
		{
			t = t + s.substr(i, 1);
			char c = nt[rd() % 4];
			t = t + c;
		}
	}

	t = t + s.substr(i - step + 1);
	return t;
}

string r_mutation(string s, double r) //Mutate with rate r
{
	int i, len = s.length();
	string t;
	for(i = 0; i < len; i++)
	{
		double z = rd() % 10000;
		if(z >= 10000 * r)
			t = t + s[i];
		else
		{
			int x = rd() % 3;

			if(x == 0)
			{
				char c = nt[(dict[s[i]] + (rd() % 3) + 1) % 4];
				t = t + c;
			}
			else if(x == 1)
			{
				t = t + s[i];
				char c = nt[rd() % 4];
				t = t + c;
			}
		}
	}

	return t;
}

int main(int argc, const char * argv[])
{
	int seqlen = stoi(argv[1]);
	double errorrate = stod(argv[2]);

	int n = atoi(argv[3]);

    dict['A'] = 0;
    dict['C'] = 1;
    dict['G'] = 2;
    dict['T'] = 3;
    
	for(int i = 0; i < n; i++)
	{
		string s = simulation(seqlen);
		string t = r_mutation(s, errorrate);
		
		printf("%s\n%s\n", s.c_str(), t.c_str());
	}
	return 0;
}
