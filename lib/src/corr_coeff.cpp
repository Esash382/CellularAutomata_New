#include "../include/corr_coeff.hpp"

uint sum(std::vector<uint> a)
{
	uint s = 0;
	for (std::size_t i = 0; i < a.size(); i++)
	{
		s += a[i];
	}
	return s;
}

double mean(std::vector<uint> a)
{
	return sum(a) / a.size();
}

uint sqsum(std::vector<uint> a)
{
	uint s = 0;
	for (std::size_t i = 0; i < a.size(); i++)
	{
		s += pow(a[i], 2);
	}
	return s;
}

double stdev(std::vector<uint> nums)
{
	uint N = nums.size();
	return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

std::vector<uint> operator-(std::vector<uint> a, uint b)
{
	std::vector<uint> retvect;
	for (std::size_t i = 0; i < a.size(); i++)
	{
		retvect.push_back(a[i] - b);
	}
	return retvect;
}

std::vector<uint> operator*(std::vector<uint> a, std::vector<uint> b)
{
	std::vector<uint> retvect;
	for (std::size_t i = 0; i < a.size() ; i++)
	{
		retvect.push_back(a[i] * b[i]);
	}
	return retvect;
}

double pearsoncoeff(std::vector<uint> X, std::vector<uint> Y)
{
	return sum((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}
