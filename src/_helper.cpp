//
// Created by Colin Small on 7/1/21.
//

#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <string>
#include <algorithm>
#include <cctype>


double factorial(double n){
    if (n == 0)
        return 1.0;
    return n * factorial(n - 1);
}


/**
 * Helper function to return multiple items from a string given a list of indices to access.
 * @param str The string to access.
 * @param indices A vector containing the indices to access and return from str.
 * @return A new vector containing only the elements at the specified indices.
 */
std::vector<char> accessMultipleIndicesString(std::string str, std::vector<int>& indices)
{
    // Make sure indices isn't larger than the array
//    assert(str.size() >= indices.size());
//    // Make sure max index in indices isn't larger than length of array
//    assert(*std::max_element(indices.begin(), indices.end()) < str.size());

    std::vector<char> returnVec = {};
    returnVec.reserve(indices.size());

    for(int i = 0; i < indices.size(); i++){
        returnVec.push_back(str[indices[i]]);
    }

    return returnVec;
}

/**
 * Returns a vector of indices where characters in a string are equal to a given character
 * @param str The string to get indices from.
 * @param comparator The character to be compared to the contents of str.
 * @return A vector of indices.
 */
std::vector<int> getIndicesWhereEqualFromString(std::string str, char comparator){
    std::vector<int> returnVec = {};

    for(int i = 0; i < str.size(); i++){
        if(str[i] == comparator)
            returnVec.push_back(i);
    }

    return returnVec;
}

double sumVector(const std::vector<double>& vec){
    double result = 0.0;
    for(auto& i : vec){
        result += i;
    }
    return result;
}

std::string toUpper(std::string s){
    for (char& c : s){
        c = toupper(c);
    }
    return s;
}

std::vector<std::string> split(const std::string& s, const std::string& delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

double meanVector(const std::vector<double>& vec){
    return sumVector(vec) / (double)vec.size();
}

// From https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring?page=1&tab=votes#tab-top
// trim from start (in place)
//static inline void ltrim(std::string &s) {
//    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
//        return !std::isspace(ch);
//    }));
//}

//// trim from end (in place)
//inline void rtrim(std::string &s) {
//    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
//        return !std::isspace(ch);
//    }).base(), s.end());
//}
//
//// trim from both ends (in place)
//inline void trim(std::string &s) {
//    ltrim(s);
//    rtrim(s);
//}
//
//// trim from start (copying)
//inline std::string ltrim_copy(std::string s) {
//    ltrim(s);
//    return s;
//}
//
//// trim from end (copying)
//inline std::string rtrim_copy(std::string s) {
//    rtrim(s);
//    return s;
//}
//
//// trim from both ends (copying)
//inline std::string trim_copy(std::string s) {
//    trim(s);
//    return s;
//}
//Ran Hu
std::string join(const std::vector<std::string> & sequence, const std::string & separator) {
    std::string result;
    for(size_t i = 0; i < sequence.size(); ++i)
        result += sequence[i] + ((i != sequence.size()-1) ? separator : "");
    return result;
}

bool cmp_pos(const std::vector<std::string> &a, const std::vector<std::string> &b) {
    return std::stoi(a[2]) < std::stoi(b[2]);
}

bool cmp_chr(const std::vector<std::string> &a, const std::vector<std::string> &b) {
    return std::stoi(a[18]) < std::stoi(b[18]);
}

bool cmp0(const std::vector<std::string> &a, const std::vector<std::string> &b) {
//    int compare = a[0].compare(b[0]);
//    if (compare == -1)
//        return true;
//    else
//        return false;
    return a[0] < b[0];
}

bool cmp1(const std::vector<std::string> &a, const std::vector<std::string> &b) {
    return std::stoi(a[1]) < std::stoi(b[1]);
}