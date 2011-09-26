/*
Copyright 2011 Frans van den Bergh. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Frans van den Bergh ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Frans van den Bergh OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of the Council for Scientific and Industrial Research (CSIR).
*/
#ifndef MTF_EDGE_QUALITY_RATING
#define MTF_EDGE_QUALITY_RATING

const double good_quality      = 1.0;
const double medium_quality    = 0.8;
const double poor_quality      = 0.5;
const double very_poor_quality = 0.2;
const double unusable_quality  = 0.0;

double edge_quality_map[46][2] = {
    {0, very_poor_quality},
    {1, poor_quality},
    {2, medium_quality},
    {3, good_quality},
    {4, good_quality},
    {5, good_quality},
    {6, good_quality},
    {7, good_quality},
    {8, good_quality},
    {9, good_quality},
    {10, good_quality},
    {11, good_quality},
    {12, medium_quality},
    {13, medium_quality},
    {14, medium_quality},
    {15, medium_quality},
    {16, medium_quality},
    {17, poor_quality},
	{18, poor_quality},    
	{19, poor_quality},
	{20, poor_quality},
	{21, poor_quality},
	{22, very_poor_quality},
	{23, very_poor_quality},
	{24, very_poor_quality},
	{25, very_poor_quality},
	{26, very_poor_quality},
	{27, poor_quality},
	{28, poor_quality},
	{29, poor_quality},
	{30, poor_quality},
	{31, poor_quality},
	{32, poor_quality},
	{33, poor_quality},
	{34, very_poor_quality},
	{35, medium_quality},
	{36, medium_quality},
	{37, medium_quality},
	{38, medium_quality},
	{39, medium_quality},
	{40, medium_quality},
	{41, medium_quality},
	{42, good_quality},
	{43, good_quality},
	{44, good_quality},
	{45, very_poor_quality},
};

#endif
