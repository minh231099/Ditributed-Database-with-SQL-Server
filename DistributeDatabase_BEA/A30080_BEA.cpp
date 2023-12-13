#include <bits/stdc++.h>

using namespace std;

ofstream fout("output.txt");

void ProcessAffinity (vector<int> &affinity, vector<vector<int>> accessFrequencyMatrix)
{
    for (int i = 0; i < accessFrequencyMatrix.size(); ++i)
    {
        int sumOfAffinity = 0;
        for (int j = 0; j < accessFrequencyMatrix[i].size(); ++j)
            sumOfAffinity += accessFrequencyMatrix[i][j];
        
        affinity[i] = sumOfAffinity;
    }
}

void ProcessAffinityMatrix (vector<vector<int>> &affinityMatrix, vector<int> affinity, vector<vector<int>> usageMatrix)
{
    for (int i = 0; i < usageMatrix.size(); ++i)
        for (int j = 0; j < usageMatrix[i].size() - 1; ++j)
            if (usageMatrix[i][j] != 0)
                for (int k = j + 1; k < usageMatrix[i].size(); ++k)
                    if (usageMatrix[i][k] != 0)
                    {
                        affinityMatrix[j][k] += affinity[i];
                        affinityMatrix[k][j] = affinityMatrix[j][k];
                    }
    // ---------------------------------------------------------------------------
    fout << "\tAffinity Matrix" << endl;
    for (int i = 0; i < affinityMatrix.size(); ++i)
        fout << "\tC" << i + 1;
    fout << endl;
    for (int i = 0; i < affinityMatrix.size(); ++i) {
        fout << "C" << i + 1 << "\t";
        for (int j = 0; j < affinityMatrix[i].size(); ++j)
            fout << affinityMatrix[i][j] << "\t";
        fout << endl;
    }
    fout << "------------------------------------------------------------" << endl;
    for (int i = 0; i < affinityMatrix.size(); ++i)
    {
        int sumOfColAffinity = 0;
        for (int j = 0; j < affinityMatrix[i].size(); ++j)
            sumOfColAffinity += affinityMatrix[i][j];
        affinityMatrix[i][i] = sumOfColAffinity;
    }
    fout << "\tAffinity Matrix sau khi đường chéo chính được tính giá trị" << endl;
    for (int i = 0; i < affinityMatrix.size(); ++i)
        fout << "\tC" << i + 1;
    fout << endl;
    for (int i = 0; i < affinityMatrix.size(); ++i) {
        fout << "C" << i + 1 << "\t";
        for (int j = 0; j < affinityMatrix[i].size(); ++j)
            fout << affinityMatrix[i][j] << "\t";
        fout << endl;
    }
    fout << "------------------------------------------------------------" << endl;
}

int Bond(vector<int> X1, vector<int> X2)
{
    int bond = 0;
    for (int i = 0; i < X1.size(); ++i)
        bond += X1[i]*X2[i];
    return bond;
}

int Cont(vector<int> X1, vector<int> X2, vector<int> X3)
{
    return 2*Bond(X1, X2) + 2*Bond(X2,X3) - 2*Bond(X1, X3);
}

void PrintPosInCont(int pos, int i, vector<vector<int> > affinityMatrix, vector<int> indexOrder) {
    indexOrder[pos] = pos;
    
    for (int j = pos; j > i; --j)
        swap(indexOrder[j], indexOrder[j-1]);

    for (int i = 0; i <= pos; ++i)
        fout << "\tC" << indexOrder[i] + 1;
    fout << endl;
    for (int i = 0; i < affinityMatrix.size(); ++i) {
        fout << "C" << i + 1;
        for (int j = 0; j <= pos; ++j)
            fout << "\t" << affinityMatrix[i][indexOrder[j]];
        fout << endl;
    }
}

void Ordering (vector<int> &newIndexOfAM, int pos, vector<vector<int>>& affinityMatrix)
{
    fout << "Tìm vị trí thích hợp của C" << pos + 1 << endl;
    int maxCont = INT_MIN;
    int newIndex = pos;
    vector<int> C0, CN;
    C0.assign(affinityMatrix.size(), 0);
    CN.assign(affinityMatrix.size(), 0);
    for (int i = 0; i <= pos; ++i) {
        int cont;
        if (i == 0)
            cont = Cont(C0, affinityMatrix[pos], affinityMatrix[newIndexOfAM[i]]);
        else if (i == pos)
            cont = Cont(affinityMatrix[newIndexOfAM[i-1]], affinityMatrix[pos], CN);
        else
            cont = Cont(affinityMatrix[newIndexOfAM[i-1]], affinityMatrix[pos], affinityMatrix[newIndexOfAM[i]]);
        PrintPosInCont(pos, i, affinityMatrix, newIndexOfAM);
        fout << "Cont: " << cont << "\n\n";
        if (maxCont < cont) {
            maxCont = cont;
            newIndex = i;
        }
    }
    fout << "Cont lớn nhất: " << maxCont << endl << "C" << pos+1 << " nằm ở vị trí " << newIndex + 1 << endl;
    newIndexOfAM[pos] = pos;
    for (int i = pos; i > newIndex; --i)
        swap(newIndexOfAM[i], newIndexOfAM[i-1]);
    fout << "----------------------------------------------------------" << endl;
}

void Reordering (vector<int> newIndex, vector<vector<int> > &affinityMatrix)
{
    vector<vector<int> > newAffinityMatrix;
    newAffinityMatrix.assign(newIndex.size(), vector<int> (newIndex.size(), 0));
    for (int i = 0; i < newIndex.size(); ++i) {
        for (int j = 0; j < newIndex.size(); ++j)
        {
            newAffinityMatrix[i][j] = affinityMatrix[newIndex[i]][newIndex[j]];
        }
    }

    affinityMatrix = newAffinityMatrix;
    fout << "\tAffinity Matrix sau khi sắp xếp lại" << endl;
    fout << "\t\t\t";
    for (int i = 0; i < affinityMatrix.size(); ++i)
        fout << "\tC" << newIndex[i] + 1;
    fout << endl;
    for (int i = 0; i < affinityMatrix.size(); ++i) {
        fout << "\t\t\tC" << newIndex[i] + 1 << "\t";
        for (int j = 0; j < affinityMatrix[i].size(); ++j)
            fout << affinityMatrix[i][j] << "\t";
        fout << endl;
    }
    fout << "------------------------------------------------------------" << endl;
}

int calculateZ (vector<int> topConner, vector<int> bottomConner, vector<vector<int> > usageMatrix, vector <int> affinity)
{   
    vector<int> ApInTCW, ApInBCW, ApInBOCW;
    int TCW = 0, BCW = 0, BOCW = 0;
    for (int i = 0; i < usageMatrix.size(); ++i) {
        vector <int> apUseC;
        bool inTop = false, inBot = false;
        for (int j = 0; j < usageMatrix[i].size(); ++j) {
            if (usageMatrix[i][j] != 0) {
                apUseC.push_back(j);
            }

            for (int k = 0; k < apUseC.size(); ++k) {
                if (topConner[apUseC[k]])
                    inTop = true;
                if (bottomConner[apUseC[k]])
                    inBot = true;
            }
        }
        if (inTop && !(inBot)) {
            TCW += affinity[i];
            ApInTCW.push_back(i);
        }
        else if (!(inTop) && inBot) { 
            BCW += affinity[i];
            ApInBCW.push_back(i);
        }
        else {
            BOCW += affinity[i];
            ApInBOCW.push_back(i);
        }
    }

    if (ApInTCW.size() == 0) {
        fout << "Không có AP nào chỉ sử dụng Top Conner" << endl;
        fout << "Nên TCW = 0" << endl;
    }
    else {
        fout << "Có ";
        for (int i = 0; i < ApInTCW.size(); i++) {
            fout << "AP" << ApInTCW[i] + 1;
            if (i != ApInTCW.size() - 1)    
                fout << ", ";
        }
        fout << " sử dụng Top Conner" << endl;
        fout << "Nên TCW = ";
        for (int i = 0; i < ApInTCW.size(); i++) {
            fout << "AFF(AP" << ApInTCW[i] + 1 << ')';
            if (i != ApInTCW.size() - 1)    
                fout << " + ";
        }
        fout << " = " << TCW << endl;
    }

    if (ApInBCW.size() == 0) {
        fout << "Không có AP nào chỉ sử dụng Bottom Conner" << endl;
        fout << "Nên BCW = 0" << endl;
    }
    else {
        fout << "Có ";
        for (int i = 0; i < ApInBCW.size(); i++) {
            fout << "AP" << ApInBCW[i] + 1;
            if (i != ApInBCW.size() - 1)    
                fout << ", ";
        }
            
        fout << " sử dụng Bottom Conner" << endl;
        fout << "Nên BCW = ";
        for (int i = 0; i < ApInBCW.size(); i++) {
            fout << "AFF(AP" << ApInBCW[i] + 1 << ')';
            if (i != ApInBCW.size() - 1)    
                fout << " + ";
        }
        fout << " = " << BCW << endl;
    }

    if (ApInBOCW.size() == 0) {
        fout << "Không có AP nào sử dụng cả Top Conner và Bottom Conner" << endl;
        fout << "Nên BOCW = 0" << endl;
    }
    else {
        fout << "Có ";
        for (int i = 0; i < ApInBOCW.size(); i++) {
            fout << "AP" << ApInBOCW[i] + 1;
            if (i != ApInBOCW.size() - 1)    
                fout << ", ";
        }
        fout << " sử dụng cả Top Conner và Bottom Conner" << endl;
        fout << "Nên BOCW = ";
        for (int i = 0; i < ApInBOCW.size(); i++) {
            fout << "AFF(AP" << ApInBOCW[i] + 1 << ')';
            if (i != ApInBOCW.size() - 1)    
                fout << " + ";
        }
        fout << " = " << BOCW << endl;
    }
    fout << "Z = " << TCW << "*" << BCW << " - " << BOCW << "^2" << " = " << TCW*BCW - BOCW*BOCW;
    return TCW*BCW - BOCW*BOCW;
}

void BEA(vector<vector<int>> &affinityMatrix, vector<vector<int>> usageMatrix, vector<int> affinity, vector<int> &left, vector<int> &right)
{
    vector<int> newIndexOfAM;
    newIndexOfAM.assign(affinityMatrix.size(), 0);
    newIndexOfAM[0] = 0;
    newIndexOfAM[1] = 1;
    for (int i = 2; i < affinityMatrix.size(); ++i)
        Ordering(newIndexOfAM, i, affinityMatrix);
    Reordering(newIndexOfAM, affinityMatrix);

    int maxZ = INT_MIN, Z;
    vector<int> topConner, bottomConner;
    topConner.assign(newIndexOfAM.size(), 1);
    bottomConner.assign(newIndexOfAM.size(), 0);
    fout << "Tính Z" << endl;
    
    for (int i = newIndexOfAM.size() - 1; i > 0; --i) {
        topConner[newIndexOfAM[i]] = 0;
        bottomConner[newIndexOfAM[i]] = 1;
        bool check = false;
        fout << "Top Conner Gồm: ";
        for (int j = 0; j < topConner.size(); ++j)
            if (topConner[newIndexOfAM[j]]) {
                if (check)
                    fout << ", ";
                fout << "C" << newIndexOfAM[j] + 1;
                check = true;
            }
        fout << endl;
        check = false;
        fout << "Bottom Conner Gồm: ";
        for (int j = 0; j < bottomConner.size(); ++j)
            if (bottomConner[newIndexOfAM[j]]) {
                if (check)
                    fout << ", ";
                fout << "C" << newIndexOfAM[j] + 1;
                check = true;
            }
        fout << endl;
        Z = calculateZ(topConner, bottomConner, usageMatrix, affinity);
        if (Z > maxZ) {
            fout << " > " << maxZ << " --> MaxZ = " << Z << endl; 
            maxZ = Z;
            left = topConner;
            right = bottomConner;
        }
        else if (Z < maxZ)
            fout << " < " << maxZ << " --> MaxZ không thay đổi" << endl;
        else
            fout << " = " << maxZ << " --> MaxZ không thay đổi" << endl;
        fout << endl;
    }

    fout << "Tính toán Inner Block" << endl;
    topConner.assign(newIndexOfAM.size(), 0);
    bottomConner.assign(newIndexOfAM.size(), 1);
    topConner[newIndexOfAM[0]] = 1;
    bottomConner[newIndexOfAM[0]] = 0;
    topConner[newIndexOfAM[newIndexOfAM.size() - 1]] = 1;
    bottomConner[newIndexOfAM[newIndexOfAM.size() - 1]] = 0;
    bool check = false;
    fout << "Top Conner Gồm: ";
    for (int j = 0; j < topConner.size(); ++j)
        if (topConner[newIndexOfAM[j]]) {
            if (check)
                fout << ", ";
            fout << "C" << newIndexOfAM[j] + 1;
            check = true;
        }
    fout << endl;
    check = false;
    fout << "Bottom Conner Gồm: ";
    for (int j = 0; j < bottomConner.size(); ++j)
        if (bottomConner[newIndexOfAM[j]]) {
            if (check)
                fout << ", ";
            fout << "C" << newIndexOfAM[j] + 1;
            check = true;
        }
    fout << endl;
    Z = calculateZ(topConner, bottomConner, usageMatrix, affinity);
    
    if (Z > maxZ) {
        fout << " > " << maxZ << " --> MaxZ = " << Z << endl; 
        maxZ = Z;
        left = topConner;
        right = bottomConner;
    }
    else if (Z < maxZ)
        fout << " < " << maxZ << " --> MaxZ không thay đổi" << endl;
    else
        fout << " = " << maxZ << " --> MaxZ không thay đổi" << endl;
    fout << "------------------------------------------------------------" << endl;
    fout << "Giá trị lớn nhất của Z là " << maxZ << " nên ta phân thành hai vùng: ";
}

int main() {
    fstream fin("test.txt");
    int ap, c, s;
    fin >> ap;
    fin >> c;
    fin >> s;

    vector<vector<int>> usageMatrix;
    usageMatrix.assign(ap, vector<int>(c, 0));

    vector<vector<int>> accessFrequencyMatrix;
    accessFrequencyMatrix.assign(ap, vector<int>(s, 0));

    for (int i = 0; i < ap; ++i)
        for (int j = 0; j < c; ++j)
            fin >> usageMatrix[i][j];
    for (int i = 0; i < ap; ++i)
        for (int j = 0; j < s; ++j)
            fin >> accessFrequencyMatrix[i][j];
    fin.close();
    
    cout << "So luong ap: " << ap << endl;
    cout << "So luong cot: " << c << endl;
    cout << "So luong may tram(s): " << s << endl;
    
    cout << "\tMa Tran Tan Suat Su Dung" << endl;
    for (int i = 0; i < ap; ++i) {
        for (int j = 0; j < c; ++j)
            cout << usageMatrix[i][j] << "\t";
        cout << endl;
    }
    cout << "\t Ma Tran Tan Suat Truy Suat" << endl;
    for (int i = 0; i < ap; ++i) {
        for (int j = 0; j < s; ++j)
            cout << accessFrequencyMatrix[i][j] << "\t";
        cout << endl;
    }
    vector<int> affinity;
    affinity.assign(ap, 0);
    ProcessAffinity(affinity, accessFrequencyMatrix);

    for (int i = 0; i < c; ++i)
        fout << "\tC" << i + 1;
    fout << "\tAffinity" << endl;
    for (int i = 0; i < ap; ++i) {
        fout << "AP" << i + 1 << '\t';
        for (int j = 0; j < c; ++j)
            fout << usageMatrix[i][j] << "\t";
        fout << affinity[i] << endl;
    }
    fout << endl;
    vector<vector<int>> affinityMatrix;
    affinityMatrix.assign(c, vector<int>(c, 0));

    ProcessAffinityMatrix(affinityMatrix, affinity, usageMatrix);
    
    vector<int> left, right;
    BEA(affinityMatrix, usageMatrix, affinity, left, right);

    fout << "(";
    bool checkWritten = false;
    for (int i = 0; i < left.size(); ++i)
    {   
        if (left[i] != 0) {
            if (checkWritten)
                fout << ", ";
            fout << "C" << i + 1;
            checkWritten = true;
        }
    }
    fout << ")";

    checkWritten = false;
    fout << " va (";
    for (int i = 0; i < right.size(); ++i)
    {   
        if (right[i] != 0) {
            if (checkWritten)
                fout << ", ";
            fout << "C" << i + 1;
            checkWritten = true;
        }
    }
    fout << ")" << endl;

    fout.close();

    return 0;
}