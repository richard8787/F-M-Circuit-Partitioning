#include "partitioner.h"
#include "cell.h"
#include "net.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;

void Partitioner::parseInput(fstream &inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str)
    {
        if (str == "NET")
        {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName)
            {
                if (cellName == ";")
                {
                    tmpCellName = "";
                    break;
                }
                else
                {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0)
                    {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else
                    {
                        if (cellName != tmpCellName)
                        {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}

void Partitioner::logicAffinity()
{
    bool allUnlock = true;
    bool partToggle = false;

    // stage 1: logic affinity
    for (int i = 0; i < _netArray.size(); i++)
    {
        vector<int> cl = _netArray[i]->getCellList();
        allUnlock = true;
        for (int j = 0; j < cl.size(); j++)
        {
            if (_cellArray[cl[j]]->getLock())
            {
                allUnlock = false;
                break;
            }
        }
        if (allUnlock)
        {
            for (int j = 0; j < cl.size(); j++)
            {
                _cellArray[cl[j]]->setPart(partToggle);
                _cellArray[cl[j]]->lock();
            }
            partToggle = !partToggle;
        }
    }

    countPartsize();

    // stage 2: move back with not meet FM balance, and force them to be half and half
    if ((!Abalance()) && (!Bbalance()))
    {
        bool more, less;
        int gap; // the gap need to move from more to less
        bool meetGap = false;
        if (getPartSize(0) > getPartSize(1))
        {
            gap = (getPartSize(0) - getPartSize(1)) / 2;
            more = 0;
            less = 1;
        }
        else
        {
            gap = (getPartSize(1) - getPartSize(0)) / 2;
            more = 1;
            less = 0;
        }

        for (int i = _netArray.size() - 1; i >= 0; i--)
        {
            vector<int> cl = _netArray[i]->getCellList();
            for (int j = 0; j < cl.size(); j++)
            {
                if (_cellArray[cl[j]]->getPart() == more)
                {
                    _cellArray[cl[j]]->setPart(less);
                    gap--;
                }
                if (gap <= 0)
                {
                    meetGap = true;
                    break;
                }
            }
            if (meetGap)
                break;
        }
    }
}

void Partitioner::countNetPartCount()
{
    for (int i = 0; i < _netArray.size(); i++)
    {
        vector<int> cl = _netArray[i]->getCellList();
        _netArray[i]->setPartCount(0, 0);
        _netArray[i]->setPartCount(1, 0);
        for (int j = 0; j < cl.size(); j++)
        {
            _netArray[i]->incPartCount(_cellArray[cl[j]]->getPart());
        }
    }
}

void Partitioner::countCutsize()
{
    int counter = 0;
    for (int i = 0; i < _netArray.size(); i++)
    {
        if ((_netArray[i]->getPartCount(0) > 0) && (_netArray[i]->getPartCount(1) > 0))
            counter++;
    }
    _cutSize = counter;
}

void Partitioner::countPartsize()
{
    int partSizeA = 0;
    int partSizeB = 0;
    for (int i = 0; i < _cellArray.size(); i++)
    {
        if (!_cellArray[i]->getPart())
            partSizeA++;
        else
            partSizeB++;
    }
    _partSize[0] = partSizeA;
    _partSize[1] = partSizeB;
}

void Partitioner::countMaxPinNum()
{
    int maxPinNum = 0;
    for (int i = 0; i < _cellArray.size(); i++)
    {
        maxPinNum = max(maxPinNum, _cellArray[i]->getPinNum());
    }
    _maxPinNum = maxPinNum;
}

void Partitioner::countGain()
{
    // initial gain to zero
    for (int i = 0; i < _cellArray.size(); i++)
    {
        _cellArray[i]->setGain(0);
    }

    // count initial gain
    for (int i = 0; i < _netArray.size(); i++)
    {
        vector<int> cl = _netArray[i]->getCellList();

        // all in each site
        if (_netArray[i]->getPartCount(0) == 0 || _netArray[i]->getPartCount(1) == 0)
        {
            for (int j = 0; j < cl.size(); j++)
            {
                _cellArray[cl[j]]->decGain();
            }
        }

        // alone on the part
        if (_netArray[i]->getPartCount(0) == 1)
        {
            for (int j = 0; j < cl.size(); j++)
            {
                if (!_cellArray[cl[j]]->getPart())
                    _cellArray[cl[j]]->incGain();
            }
        }
        if (_netArray[i]->getPartCount(1) == 1)
        {
            for (int j = 0; j < cl.size(); j++)
            {
                if (_cellArray[cl[j]]->getPart())
                    _cellArray[cl[j]]->incGain();
            }
        }
    }
}

bool Partitioner::Abalance()
{
    double lb = (1 - getBFactor()) / 2 * (double)getCellNum();
    double ub = (1 + getBFactor()) / 2 * (double)getCellNum();
    bool meetlb = (((getPartSize(0) - 1) >= lb) && ((getPartSize(1) + 1) >= lb));
    bool meetub = (((getPartSize(0) - 1) <= ub) && ((getPartSize(1) + 1) <= ub));
    return meetlb && meetub;
}

bool Partitioner::Bbalance()
{
    double lb = (1 - getBFactor()) / 2 * (double)getCellNum();
    double ub = (1 + getBFactor()) / 2 * (double)getCellNum();
    bool meetlb = (((getPartSize(0) + 1) >= lb) && ((getPartSize(1) - 1) >= lb));
    bool meetub = (((getPartSize(0) + 1) <= ub) && ((getPartSize(1) - 1) <= ub));
    return meetlb && meetub;
}

int Partitioner::FM()
{
    int maxPartialSum = INT32_MIN;
    int maxPartialSumID = 0;
    bool choose = false;
    const int set_size = 2 * getMaxPinNum() + 1;
    vector<int> partialSum(_cellArray.size());
    vector<Cell *> moveCell(_cellArray.size());
    bool F = false, T = false;           // set FromSet and ToSet
    unordered_set<Cell *> set[set_size]; // bucket list (using unordered_set)
    Cell *move = nullptr;                // move which cell
    vector<int> nl;                      // netlist of the cell
    vector<int> cl;                      // celllist of the net

    // create initial bucket list
    for (int i = 0; i < _cellArray.size(); i++)
    {
        _cellArray[i]->unlock(); // unlock all
        set[getMaxPinNum() + _cellArray[i]->getGain()].insert(_cellArray[i]);
    }

    // move all
    for (int itt = 0; itt < _cellArray.size(); itt++)
    {
        // check can move or not
        choose = false;

        //  choose which to move
        for (int i = set_size - 1; i >= 0; i--)
        {
            if (choose)
                break;

            if (!set[i].empty())
            {
                for (auto iter = set[i].begin(); iter != set[i].end(); iter++)
                {
                    // cout << (*iter)->getName() << endl;
                    if (!(*iter)->getPart() && Abalance() && !(*iter)->getLock()) // move from partA legal
                    {
                        move = (*iter);
                        moveCell[itt] = move;
                        choose = true;
                        // cout << "move: " << move->getName() << endl;
                        break;
                    }
                    else if ((*iter)->getPart() && Bbalance() && !(*iter)->getLock()) // move from partB legal
                    {
                        move = (*iter);
                        moveCell[itt] = move;
                        choose = true;
                        // cout << "move: " << move->getName() << endl;
                        break;
                    }
                }
            }
        }

        if (!choose)
            cout << "Warning: not choose any thing!!!!!!!!!!!!!!!" << endl;

        // check the move cell's FromPart and ToPart
        F = move->getPart();
        T = !F;

        // move the cell to opposite part
        move->setPart(T);
        move->lock();
        for (int i = 0; i < set_size; i++)
            set[i].erase(move);
        _partSize[F]--;
        _partSize[T]++;

        // count the maximum partial sum and id
        if (itt == 0)
            partialSum[itt] = move->getGain();
        else
            partialSum[itt] = partialSum[itt - 1] + move->getGain();

        if (partialSum[itt] > maxPartialSum)
        {
            maxPartialSum = partialSum[itt];
            maxPartialSumID = itt;
        }

        // Update Gain
        nl = move->getNetList();
        for (int i = 0; i < nl.size(); i++)
        {
            // T
            cl = _netArray[nl[i]]->getCellList();
            if (_netArray[nl[i]]->getPartCount(T) == 0)
            {
                for (int j = 0; j < cl.size(); j++)
                {
                    if (!_cellArray[cl[j]]->getLock())
                    {
                        set[getMaxPinNum() + _cellArray[cl[j]]->getGain()].erase(_cellArray[cl[j]]);
                        _cellArray[cl[j]]->incGain();
                        set[getMaxPinNum() + _cellArray[cl[j]]->getGain()].insert(_cellArray[cl[j]]);
                    }
                }
            }
            else if (_netArray[nl[i]]->getPartCount(T) == 1)
            {
                for (int j = 0; j < cl.size(); j++)
                {
                    if (_cellArray[cl[j]]->getPart() == T)
                    {
                        if (!_cellArray[cl[j]]->getLock())
                        {
                            set[getMaxPinNum() + _cellArray[cl[j]]->getGain()].erase(_cellArray[cl[j]]);
                            _cellArray[cl[j]]->decGain();
                            set[getMaxPinNum() + _cellArray[cl[j]]->getGain()].insert(_cellArray[cl[j]]);
                        }
                    }
                }
            }

            // F(n) <- F(n)-1; T(n)<-T(n)+1;
            _netArray[nl[i]]->decPartCount(F);
            _netArray[nl[i]]->incPartCount(T);

            // F
            if (_netArray[nl[i]]->getPartCount(F) == 0)
            {
                for (int j = 0; j < cl.size(); j++)
                {
                    if (!_cellArray[cl[j]]->getLock())
                    {
                        set[getMaxPinNum() + _cellArray[cl[j]]->getGain()].erase(_cellArray[cl[j]]);
                        _cellArray[cl[j]]->decGain();
                        set[getMaxPinNum() + _cellArray[cl[j]]->getGain()].insert(_cellArray[cl[j]]);
                    }
                }
            }
            else if (_netArray[nl[i]]->getPartCount(F) == 1)
            {
                for (int j = 0; j < cl.size(); j++)
                {
                    if (_cellArray[cl[j]]->getPart() == F)
                    {
                        if (!_cellArray[cl[j]]->getLock())
                        {
                            set[getMaxPinNum() + _cellArray[cl[j]]->getGain()].erase(_cellArray[cl[j]]);
                            _cellArray[cl[j]]->incGain();
                            set[getMaxPinNum() + _cellArray[cl[j]]->getGain()].insert(_cellArray[cl[j]]);
                        }
                    }
                }
            }
        }

        // show the set
        // for (int i = 0; i < set_size; i++)
        // {
        //     cout << "gain =" << i - getMaxPinNum() << " : ";
        //     for (const auto &s : set[i])
        //     {
        //         cout << s->getName() << " ";
        //     }
        //     cout << endl;
        // }
        // reportCellGain();
        // reportCellPart();
        // cout << "-----------------------------move a cell end-----------------------------" << endl;
    }

    // move back
    if (maxPartialSum > 0)
    {
        for (int i = maxPartialSumID + 1; i < moveCell.size(); i++)
        {
            moveCell[i]->setPart(!moveCell[i]->getPart());
        }
    }

    countNetPartCount();
    countPartsize();
    countCutsize();
    countGain();

    // show
    // for (int i = 0; i < _cellArray.size(); i++)
    // {
    //     cout << "partialSum: " << partialSum[i] << " move: " << moveCell[i]->getName() << endl;
    // }
    // cout << "maxPartialSum: " << maxPartialSum << endl;
    // cout << "-----------------------------one iteration end-----------------------------" << endl;

    return maxPartialSum;
}

void Partitioner::partition()
{
    // set the initial partition
    logicAffinity();
    countNetPartCount();
    countPartsize();
    countCutsize();
    countGain();
    countMaxPinNum();

    // report initial partition
    cout << "****initial partition****" << endl;
    reportMaxPinNum();
    reportCutsize();
    cout << endl;
    // reportNetPartCount();
    // reportCellPart();
    // reportCellGain();

    // do the FM
    int MPS = 777; // maximum partial sum in each itheration
    const int earlyBreakTime = 100;
    vector<int> count5(5);
    for (int i = 1; (i <= 150) && (MPS > 0); i++)
    {
        cout << "****iteration: " << i << "****" << endl;
        MPS = FM();
        cout << "maxPartialSum: " << MPS << endl;
        reportCutsize();
        cout << endl;

        count5[i % 5] = MPS;
        if ((i > earlyBreakTime) && (accumulate(count5.begin(), count5.end(), 0) < 20))
        {
            cout << "Early Break" << endl;
            break;
        }
    }

    // final report
    // reportNetPartCount();
    // reportCellPart();
    // reportCellGain();
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i)
    {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j)
        {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i)
    {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j)
        {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCellGain() const
{
    cout << "==================== Cell Gain ====================" << endl;
    for (int i = 0, end_i = _cellArray.size(); i < end_i; ++i)
    {
        cout << setw(8) << _cellArray[i]->getName() << " gain: " << _cellArray[i]->getGain() << endl;
    }
}

void Partitioner::reportCellPart() const
{
    cout << "==================== Cell Part ====================" << endl;
    for (int i = 0, end_i = _cellArray.size(); i < end_i; ++i)
    {
        cout << setw(8) << _cellArray[i]->getName() << " part: " << _cellArray[i]->getPart() << endl;
    }
}

void Partitioner::reportNetPartCount() const
{
    cout << "==================== Net PartCount ====================" << endl;
    for (int i = 0; i < _netArray.size(); i++)
    {
        cout << _netArray[i]->getName() << ": " << _netArray[i]->getPartCount(0) << " " << _netArray[i]->getPartCount(1) << endl;
    }
}

void Partitioner::reportCutsize() const
{
    cout << "CutSize is: " << getCutSize() << endl;
}

void Partitioner::reportMaxPinNum() const
{
    cout << "maxPinNum is: " << getMaxPinNum() << endl;
}

void Partitioner::writeResult(fstream &outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i)
    {
        if (_cellArray[i]->getPart() == 0)
        {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i)
    {
        if (_cellArray[i]->getPart() == 1)
        {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i)
    {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i)
    {
        delete _netArray[i];
    }
    return;
}