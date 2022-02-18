
clear all, clc

data=xlsread('DATA.xlsx');

results=NbFactors(data,20);