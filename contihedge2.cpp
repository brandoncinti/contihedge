// contihedge3.cpp : Defines the entry point for the console application.
// This funciton is a single thread implementation of continuous hedging paper by Ederman.
// Author: Xin Wang
// 
#include "stdafx.h"
#include <cmath>
#include <math.h>
#include <boost\math\distributions.hpp>
#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include "contihedge3.hpp"

using namespace std;

ContHedge::ContHedge() : nsims(10000), intr(0.05), vsigma(0.20), ytime(0.083333), 
	c0(2.512), s0(100.00), sk(100.00)
{
	// initialize the given constants in constructor;
	cout << "Initialization success" << endl;

}

ContHedge::ContHedge(const ContHedge& che): nsims(che.nsims), intr(che.intr), 
	vsigma(che.vsigma), ytime(che.ytime), c0(che.c0), s0(che.s0), sk(che.sk)
{
	// Copy constructor
}

ContHedge::~ContHedge(){
	cout << "Finish the hedge simulation." << endl;
}

ContHedge& ContHedge::operator = (const ContHedge& camt)
{
	// Assignment operator definition;
	if (this == &camt)
		return *this;
}

double ContHedge::optiondelta(double tau, double si)
{
	double d1 = (std::log(si/sk)+(intr+0.5*vsigma*vsigma)*tau)/(vsigma*std::sqrt(tau));
	double d2 = d1 - vsigma*std::sqrt(tau); 
	boost::math::normal_distribution<> Ncdf(0,1);
	double k = 1.0/(1.0 + 0.2316419*d1);
	double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
	if (TypeOption == 'C')
		return cdf(Ncdf, d1);
		// return (1.0 - (1.0/(pow(2*M_PI, 0.5)))*exp(-0.5*d1*d1) * k_sum);
	else
		return (1.0 - (1.0/(pow(2*M_PI, 0.5)))*exp(-0.5*d1*d1) * k_sum) - 1;
}

double ContHedge::optionvalue(double tau, double si)
{
	double d1 = (std::log(si/sk)+(intr+0.5*vsigma*vsigma)*tau)/(vsigma*std::sqrt(tau));
	double d2 = d1 - vsigma*std::sqrt(tau); 
	boost::math::normal_distribution<> Ncdf(0,1);
	if (TypeOption == 'C')
		return si*cdf(Ncdf, d1)-sk*exp(-1*intr*tau)*cdf(Ncdf, d2);// call option value;
	else
		return sk*exp(-1*intr*tau)*(1-cdf(Ncdf, d2)) - si*(1-cdf(Ncdf, d1)); // put option value;
}

void ContHedge::userinit()
{
	// Prompt input function to initialize parameters;
	cout << "Please input hedging paratmeters: \n"<< endl;
	cout << "Number of simulations: " << std::endl;
	cin >> nsims;

	cout << "Interest Rate: " << endl;
	cin >> intr;

	cout << "Stock Appreciate Rate: " << endl;
	cin >> stckr;

	cout << "Volatility Sigma: " << endl;
	cin >> vsigma;

	cout << "Mature time: " << endl;
	cin >> ytime;

	cout << "Intial Option price: " << endl;
	cin >> c0;

	cout << "Strike price: " << endl;
	cin >> sk;

	cout << "Option Type: " << endl;
	cin >> TypeOption;

	cout << "Initial stock price: " << endl;
	cin >> s0;
	
	cout << "Parameters replaced by user input values" << endl;

}

void ContHedge::deltahedge(int numtrades, std::vector<double>& finalporl) 
{
	// delta hedge for both put and call
	// stock/option = delta_call;
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0, vsigma); // mean 0, sigma volatile;		
	// individual stock price simulation;
	for (int indsim=0; indsim<nsims; indsim++)
	{		
		// Part 1: initialize trade with known call option price
		double cashbalance;
		double tau = ytime;
		double nsteps = numtrades;
		double dt = ytime/nsteps;
		double p0 = sk*exp(-1*intr*tau)+c0-s0;
		cashbalance = optionvalue(tau, s0);
		// cout << "The inital option value is :" << cashbalance;		
		double stockprice = s0; 		
		double stockhold = optiondelta(tau, stockprice); 		
		cashbalance -= stockhold*s0;
		
		// Part 2: Re-hedging with discrete intevals;
		for (int hstep = 1; hstep<nsteps; hstep++){
			// generate the stock price; 
			double tau = (nsteps-hstep)*dt; 					
			double nnumrand = distribution(generator); 
			stockprice *= exp((stckr-0.5*vsigma*vsigma)*dt+vsigma*sqrt(dt)*nnumrand);
			// stockprice = s0*exp(vsigma*sqrt(hstep*dt)*nnumrand + (stckr-0.5*vsigma*vsigma)*hstep*dt);		
			// rehedge: delta;
			double hdelta = optiondelta(tau, stockprice); 
			cashbalance *= exp(stckr*dt);			
			cashbalance -= (hdelta - stockhold)*stockprice;
			stockhold = hdelta; // new stock holding amount;
			// std::cout << "time to end is " << tau << std::endl;
			// std::cout << "stockprice is " << stockprice << std::endl;
		}
		distribution.reset();

		// Part 3: finalize trade;
		// stock final value;	
		double NetPosition = 0;
		double nnumrand = distribution(generator);		
		// stockprice *= 1+stckr*dt+vsigma*sqrt(dt)*nnumrand; 
		// stockprice *= exp((stckr-0.5*vsigma*vsigma)*dt+vsigma*sqrt(dt)*nnumrand);
		stockprice = s0*exp(vsigma*sqrt(ytime)*nnumrand + (stckr-0.5*vsigma*vsigma)*ytime);		
		cashbalance *= exp(stckr*dt);				
		// double hdelta = optiondelta(0, stockprice); 
		// cashbalance -= (hdelta-stockhold)*stockprice; // last delta: time left 0;
		if (TypeOption == 'C'){
			// Sell a Call option 	
			double EOpayoff = stockprice > sk ? (stockprice - sk) : 0.0;		
			// double EOpayoff = stockprice > sk ? sk : 0.0;
			NetPosition = cashbalance - EOpayoff + stockhold*stockprice ; // pay money, payoff, sell stock;
			// NetPosition = cashbalance + EOpayoff ; // 
		}
		else		
		{
			// Sell a Put option
			// double EOpayoff = stockprice < sk ? sk : 0.0;		
			// double final_stockhold = stockprice < sk ? -1 : 0.0;
			// cashbalance -= (final_stockhold+stockhold)*stockprice; // cash spending update;
			double EOpayoff = stockprice < sk ? (sk-stockprice) : 0.0;
			NetPosition = cashbalance - EOpayoff + stockhold*stockprice;
		}
		// Return money ammount;
		finalporl[indsim] = NetPosition;

	}
	distribution.reset();
	cout << "Finished deltahedge function" << std::endl;
}

