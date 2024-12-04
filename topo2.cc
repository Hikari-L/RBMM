#include "common.h"
#include "send.h"

/*
struct FlowInfo {
    std::string eid;
    int pktCount;
    int numPaths;
    int pktSize;
};

std::vector<FlowInfo> ReadFlowDataFromFile(const std::string &filename) {
    std::vector<FlowInfo> flowData;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return flowData;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        FlowInfo flow;
        std::string temp;

        if (std::getline(iss, temp, ' ')) flow.eid = temp.substr(5); 
        std::getline(iss, temp, ' ');                                
        std::getline(iss, temp, ' '); flow.pktCount = std::stoi(temp);
        std::getline(iss, temp, ' ');                              
        std::getline(iss, temp, ' ');                             
        std::getline(iss, temp, ' '); flow.numPaths = std::stoi(temp);
        std::getline(iss, temp, ' ');                               
        std::getline(iss, temp, ' ');                               
        flow.pktSize = std::stoi(temp);                             

        flowData.push_back(flow);
    }

    file.close();
    return flowData;
}

*/


struct RelayPath
{
    Address addr;
    string bandwidthMbps;
};

std::string stringAddition(const std::string &a, const std::string &b)
{
    int i = a.size() - 1, j = b.size() - 1, carry = 0;
    std::string result = "";
    while (i >= 0 || j >= 0 || carry)
    {
        int sum = carry;
        if (i >= 0)
            sum += a[i--] - '0';
        if (j >= 0)
            sum += b[j--] - '0';
        result += (sum % 10) + '0';
        carry = sum / 10;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

uint32_t selectNum(const std::vector<uint32_t> &counts, const std::vector<double> &weights)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> dist(weights.begin(), weights.end());
    return counts[dist(gen)];
}

uint32_t poissonRandom(double mean)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<> dist(mean);
    return dist(gen);
}

struct PacketCompare
{
    bool operator()(const Ptr<Packet> &lhs, const Ptr<Packet> &rhs) const
    {
        IDPHeader leftidp;
        lhs->RemoveHeader(leftidp);
        uint32_t l_seq = leftidp.GetSeqNum();
        IDPHeader rightidp;
        rhs->RemoveHeader(rightidp);
        uint32_t r_seq = rightidp.GetSeqNum();
        lhs->AddHeader(leftidp);
        rhs->AddHeader(rightidp);
        return l_seq < r_seq;
    }
};

struct ScoreEntry
{
    double score;
    string srcEid;
    bool operator<(const ScoreEntry &other) const
    {
        return score > other.score || (score == other.score && srcEid > other.srcEid);
    }
};

set<ScoreEntry> scoreSet;

void PrintScoreSet(const std::set<ScoreEntry> &scoreSet)
{
    std::cout << "Current scores in scoreSet:" << std::endl;
    for (const auto &entry : scoreSet)
    {
        std::cout << "Eid: " << entry.srcEid << ", Score: " << entry.score << std::endl;
    }
}

void updateScore(const string &srcEid, double deltaScore)
{
    bool found = false;
    for (auto it = scoreSet.begin(); it != scoreSet.end(); ++it)
    {
        if (it->srcEid == srcEid)
        {
            double score = it->score;
            scoreSet.erase(it);
            scoreSet.insert({score + deltaScore, srcEid});
            found = true;
            break;
        }
    }
    if (!found)
    {
        scoreSet.insert({deltaScore, srcEid});
    }
}

void Bzeroscore(const string &srcEid, double Score)
{
    for (auto it = scoreSet.begin(); it != scoreSet.end(); ++it)
    {
        if (it->srcEid == srcEid)
        {
            scoreSet.erase(it);
            break;
        }
    }
    scoreSet.insert({Score, srcEid});
}

string getMaxScoreSrcEid()
{
    if (!scoreSet.empty())
    {
        return scoreSet.begin()->srcEid;
    }
    return "";
}

string getMinScoreSrcEid()
{
    if (!scoreSet.empty())
    {
        return scoreSet.rbegin()->srcEid;
    }
    return "";
}

double getScoreBySrcEid(const string &srcEid)
{
    for (const auto &entry : scoreSet)
    {
        if (entry.srcEid == srcEid)
        {
            return entry.score;
        }
    }
    return -1;
}

void calScore(uint32_t SeqNum, const std::string &srcEid,
              std::map<std::string, std::set<Ptr<Packet>, PacketCompare>> &bufferQueues,
              std::map<std::string, uint32_t> &expectedSeqNums, double Lnorm, int flag1)
{
    if (bufferQueues.find(srcEid) != bufferQueues.end() && !bufferQueues[srcEid].empty())
    {
        Ptr<Packet> firstPacket = *bufferQueues[srcEid].begin();
        IDPHeader firstIdp;
        firstPacket->RemoveHeader(firstIdp);
        uint32_t firstSeq = firstIdp.GetSeqNum();
        firstPacket->AddHeader(firstIdp);

        auto lastElementIterator = --bufferQueues[srcEid].end();
        Ptr<Packet> lastpacket = *lastElementIterator;
        IDPHeader lastIdp;
        lastpacket->RemoveHeader(lastIdp);
        uint32_t lastSeq = lastIdp.GetSeqNum();
        lastpacket->AddHeader(lastIdp);

        uint32_t D1 = firstSeq - expectedSeqNums[srcEid];
        uint32_t Dinside = lastSeq - firstSeq - bufferQueues[srcEid].size() + 1;

        if (flag1 == 0)
        {
            if (SeqNum >= lastSeq + 1)
            {
                updateScore(srcEid, static_cast<double>(D1 + Dinside) * Lnorm);
            }
            else
            {
                updateScore(srcEid, static_cast<double>(D1 + Dinside - lastSeq + SeqNum - 1) * Lnorm);
            }
        }
        else if (flag1 == 1)
        {
            updateScore(srcEid, static_cast<double>(-1 * D1) * Lnorm);
        }
    }
    else
    {
        cout << "BufferQueue for " << srcEid << " is empty or does not exist." << endl;
        scoreSet.insert({0, srcEid});
    }
}

double generateParetoFlowSize(double alpha, double xm, double maxFlowSize)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double u = dis(gen);
    double flowSize = xm / std::pow(1 - u, 1.0 / alpha);

    if (flowSize > maxFlowSize)
    {
        flowSize = maxFlowSize;
    }

    return flowSize;
}


double CalculateProbabilityForEid(const std::string &eid, const std::set<ScoreEntry> &scoreSet)
{
    double minLogScore = std::numeric_limits<double>::infinity();
    double maxLogScore = -std::numeric_limits<double>::infinity();
    double curLogScore = 0;
    bool found = false;

    for (const auto &scoreEntry : scoreSet)
    {
        double logScore = std::log(scoreEntry.score + 1);
        if (logScore < minLogScore)
        {
            minLogScore = logScore;
        }
        if (logScore > maxLogScore)
        {
            maxLogScore = logScore;
        }
        if (scoreEntry.srcEid == eid)
        {
            curLogScore = logScore;
            found = true;
        }
    }

    if (!found)
    {
        std::cerr << "Error: Eid not found in scoreSet." << std::endl;
        return -1.0;
    }
    if (maxLogScore == minLogScore)
    {
        return 0.0;
    }

    double epsilon = 1e-5;
    if ((maxLogScore - minLogScore > std::log(1e6)) && std::abs(curLogScore - minLogScore) < epsilon)
    {
        return 1.0;
    }

    double probability = (curLogScore - minLogScore) / (maxLogScore - minLogScore);
    return probability;
}

int IsPersistLargeFlow(const std::string &srcEid, uint32_t Fsizecur, uint32_t Psize, uint32_t CurStandbyBuffer, std::map<std::string, bool> &Mk,
                       std::map<std::string, Time> &Rk, Time now)
{
    auto mk = Mk.find(srcEid);
    auto rk = Rk.find(srcEid);

    if ((Fsizecur * Psize > alpha * CurStandbyBuffer) && (mk != Mk.end() && mk->second == true))
    {
        if (rk != Rk.end())
        {
            Time timeDifference = now - rk->second;
            if (timeDifference.GetMilliSeconds() > T && timeDifference.GetMilliSeconds() < 100)
            {
                return 1;
            }
            else if (timeDifference.GetMilliSeconds() >= 100)
            {
                return 2;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            cout << "Err: Can't find rk for this Eid!" << endl;
            return 0;
        }
    }
    else if ((Fsizecur * Psize > alpha * CurStandbyBuffer) && (mk != Mk.end() && mk->second == false))
    {
        Mk[srcEid] = true;
        if (rk != Rk.end())
        {
            Rk[srcEid] = now;
            return 0;
        }
        else
        {
            cout << "Err: Can't find rk for this Eid!" << endl;
            return 0;
        }
    }
    else
    {
        return 0;
    }
}

void PrintBufferQueue(std::map<std::string, std::set<Ptr<Packet>, PacketCompare>> bufferQueues)
{
    for (const auto &entry : bufferQueues)
    {
        const std::string &eid = entry.first;
        const std::set<Ptr<Packet>, PacketCompare> &packetSet = entry.second;
        std::cout << "EID: " << eid << " has " << packetSet.size() << " packets remaining." << std::endl;

        for (const auto &packet : packetSet)
        {
            std::cout << "Packet with size: " << packet->GetSize() << std::endl;
        }
    }
}

void transit(Ptr<Socket> sock)
{
    Ptr<Packet> packet = sock->Recv();
    Ipv6Header ipHdr;
    packet->RemoveHeader(ipHdr);
    IDPHeader idp;
    packet->RemoveHeader(idp);
    std::string dstEid = idp.GetDestinationEid().substr(0, 20);
    std::string srcEid = idp.GetSourceEid().substr(0, 20);
    uint32_t seq = idp.GetSeqNum();
    uint8_t type = idp.GetType();
    std::string dstEid = idp.GetDestinationEid().substr(0, 20);
    Ipv6Address dstAddr = eidInquire(dstEid)[0];

    if (type == 'R')
    {
        Ptr<Node> n = sock->GetNode();
        Ptr<Socket> retransSock = Socket::CreateSocket(n, Ipv6RawSocketFactory::GetTypeId());
        retransSock->Bind();
        retransSock->Connect(Inet6SocketAddress(dstAddr));
        retransSock->SetAttribute("Protocol", UintegerValue(153));
        packet->AddHeader(idp);
        retransSock->Send(packet);
        retransSock->Close();
    }
    else if (type == 'D')
    {
        packet->AddHeader(idp);
        sock->SendTo(packet, 0, Inet6SocketAddress(dstAddr));
    }
    else
    {
        packet->AddHeader(idp);
        sock->SendTo(packet, 0, Inet6SocketAddress(dstAddr));
    }
}

void CS(Ptr<Socket> sock)
{
    static uint32_t globalBufferQueue;
    static std::map<std::string, std::set<Ptr<Packet>, PacketCompare>> bufferQueues;
    static std::map<std::string, uint32_t> expectedSeqNums;
    Address boundAddress;
    sock->GetSockName(boundAddress);
    Inet6SocketAddress reorderAddr = Inet6SocketAddress::ConvertFrom(boundAddress);
    if (reorderAddr.GetIpv6() == Ipv6Address("2001:3::200:ff:fe00:5"))
    {
        Ptr<Packet> packet = sock->Recv();
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);
        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();
        uint8_t type = idp.GetType();
        static std::map<std::string, ns3::Time> eidArrivalTimes;
        if (expectedSeqNums.find(srcEid) == expectedSeqNums.end())
        {
            expectedSeqNums[srcEid] = 1;
        }

        if (seq == expectedSeqNums[srcEid])
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
            expectedSeqNums[srcEid]++;
            while ((bufferQueues.find(srcEid) != bufferQueues.end()) && (!bufferQueues[srcEid].empty()))
            {
                auto pktIterator = bufferQueues[srcEid].begin();
                Ptr<Packet> bufferedPacket = *pktIterator;
                IDPHeader bufferedIdp;
                bufferedPacket->RemoveHeader(bufferedIdp);
                if (bufferedIdp.GetSeqNum() == expectedSeqNums[srcEid])
                {
                    uint32_t pktSize;
                    pktSize = bufferedIdp.GetPacketSize();
                    bufferedPacket->AddHeader(bufferedIdp);
                    sock->SendTo(bufferedPacket, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
                    expectedSeqNums[srcEid]++;
                    globalBufferQueue -= pktSize;
                    bufferQueues[srcEid].erase(pktIterator);

                    if (bufferQueues[srcEid].empty())
                    {
                        bufferQueues.erase(srcEid);
                    }
                }
                else
                {
                    bufferedPacket->AddHeader(bufferedIdp);
                    break;
                }
            }
        }

        else if (seq > expectedSeqNums[srcEid])
        {
            uint32_t Psize = idp.GetPacketSize();
            if (bufferQueues.find(srcEid) == bufferQueues.end())
            {
                bufferQueues[srcEid] = std::set<Ptr<Packet>, PacketCompare>();
                eidArrivalTimes[srcEid] = ns3::Simulator::Now();
            }
            if (maxBufferSize <= globalBufferQueue)
            {
                cout << "full." << endl;
                cout << "Eid : " << srcEid << " 's " << idp.GetSeqNum() << "packet." << endl;
                packet->AddHeader(idp);
                int sendbyte = sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
                cout << "send: " << sendbyte << endl;
            }
            else
            {
                packet->AddHeader(idp);
                auto res = bufferQueues[srcEid].insert(packet);
                if (res.second)
                {
                    globalBufferQueue += Psize;
                }
            }
        }
        else
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
        }
    }
    else if (reorderAddr.GetIpv6() == Ipv6Address("2001:4::200:ff:fe00:7"))
    {
        Ptr<Packet> packet = sock->Recv();
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);
        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();
        uint8_t type = idp.GetType();
        static std::map<std::string, ns3::Time> eidArrivalTimes;
        if (expectedSeqNums.find(srcEid) == expectedSeqNums.end())
        {
            expectedSeqNums[srcEid] = 1;
        }
        if (seq == expectedSeqNums[srcEid])
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:8"));

            expectedSeqNums[srcEid]++;
            while ((bufferQueues.find(srcEid) != bufferQueues.end()) && (!bufferQueues[srcEid].empty()))
            {
                auto pktIterator = bufferQueues[srcEid].begin();
                Ptr<Packet> bufferedPacket = *pktIterator;
                IDPHeader bufferedIdp;
                bufferedPacket->RemoveHeader(bufferedIdp);
                if (bufferedIdp.GetSeqNum() == expectedSeqNums[srcEid])
                {
                    uint32_t pktSize;
                    pktSize = bufferedIdp.GetPacketSize();
                    bufferedPacket->AddHeader(bufferedIdp);
                    sock->SendTo(bufferedPacket, 0, Inet6SocketAddress("2001:3::200:ff:fe00:8"));
                    expectedSeqNums[srcEid]++;
                    globalBufferQueue -= pktSize;
                    bufferQueues[srcEid].erase(pktIterator);
                    if (bufferQueues[srcEid].empty())
                    {
                        bufferQueues.erase(srcEid);
                    }
                }
                else
                {
                    bufferedPacket->AddHeader(bufferedIdp);
                    break;
                }
            }
        }

        else if (seq > expectedSeqNums[srcEid])
        {
            uint32_t Psize = idp.GetPacketSize();
            if (bufferQueues.find(srcEid) == bufferQueues.end())
            {
                bufferQueues[srcEid] = std::set<Ptr<Packet>, PacketCompare>();
                eidArrivalTimes[srcEid] = ns3::Simulator::Now();
            }
            if (maxBufferSize <= globalBufferQueue)
            {
                cout << "full." << endl;
                cout << "Eid : " << srcEid << " 's " << idp.GetSeqNum() << "packet." << endl;
                packet->AddHeader(idp);
                int sendbyte = sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:8"));
                cout << "send: " << sendbyte << endl;
            }
            else
            {
                packet->AddHeader(idp);
                auto res = bufferQueues[srcEid].insert(packet);
                if (res.second)
                {
                    globalBufferQueue += Psize;
                }
            }
        }
        else
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:8"));
        }
    }
}

void NOREORDER(Ptr<Socket> sock)
{
    Address boundAddress;
    sock->GetSockName(boundAddress);
    Inet6SocketAddress reorderAddr = Inet6SocketAddress::ConvertFrom(boundAddress);
    if (reorderAddr.GetIpv6() == Ipv6Address("2001:3::200:ff:fe00:5"))
    {
        Ptr<Packet> packet = sock->Recv();
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);
        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();
        packet->AddHeader(idp);
        sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
    }
    else if (reorderAddr.GetIpv6() == Ipv6Address("2001:4::200:ff:fe00:7"))
    {
        Ptr<Packet> packet = sock->Recv();
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);
        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();
        packet->AddHeader(idp);
        sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:8"));
    }
}

void RCBRB(Ptr<Socket> sock)
{

    static uint32_t globalBufferQueue1 = 0;
    static uint32_t globalBufferQueue2 = 0;
    static uint32_t maxGlobalBufferQueue1 = 0;
    static uint32_t maxGlobalBufferQueue2 = 0;
    static std::map<std::string, std::set<Ptr<Packet>, PacketCompare>> bufferQueues;
    static std::map<std::string, uint32_t> expectedSeqNums;
    Address boundAddress;
    sock->GetSockName(boundAddress);
    Inet6SocketAddress reorderAddr = Inet6SocketAddress::ConvertFrom(boundAddress);
    Time firstReceiveTime = Simulator::Now();
    static Time lastReceiveTime = Seconds(0);
    if (firstReceiveTime > lastReceiveTime + Seconds(0.5))
    {
        outfile4 << "QUEUE1: " << maxGlobalBufferQueue1 << "\n";
        outfile4 << "QUEUE2: " << maxGlobalBufferQueue2 << "\n";
        lastReceiveTime = firstReceiveTime;
    }
    if (reorderAddr.GetIpv6() == Ipv6Address("2001:3::200:ff:fe00:5"))
    {
        Ptr<Packet> packet = sock->Recv();
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);
        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();
        uint8_t type = idp.GetType();
        static std::map<std::string, ns3::Time> eidArrivalTimes;
        if (expectedSeqNums.find(srcEid) == expectedSeqNums.end())
        {
            expectedSeqNums[srcEid] = 1;
        }
        if (seq == expectedSeqNums[srcEid])
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
            expectedSeqNums[srcEid]++;
            while ((bufferQueues.find(srcEid) != bufferQueues.end()) && (!bufferQueues[srcEid].empty()))
            {
                auto pktIterator = bufferQueues[srcEid].begin();
                Ptr<Packet> bufferedPacket = *pktIterator;
                IDPHeader bufferedIdp;
                bufferedPacket->RemoveHeader(bufferedIdp);
                if (bufferedIdp.GetSeqNum() == expectedSeqNums[srcEid])
                {
                    uint32_t pktSize;
                    pktSize = bufferedIdp.GetPacketSize();
                    bufferedPacket->AddHeader(bufferedIdp);
                    sock->SendTo(bufferedPacket, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
                    expectedSeqNums[srcEid]++;
                    globalBufferQueue1 -= pktSize;
                    bufferQueues[srcEid].erase(pktIterator);

                    if (bufferQueues[srcEid].empty())
                    {
                        bufferQueues.erase(srcEid);
                    }
                }
                else
                {
                    bufferedPacket->AddHeader(bufferedIdp);
                    break;
                }
            }
        }

        else if (seq > expectedSeqNums[srcEid])
        {
            uint32_t Psize = idp.GetPacketSize();
            if (bufferQueues.find(srcEid) == bufferQueues.end())
            {
                bufferQueues[srcEid] = std::set<Ptr<Packet>, PacketCompare>();
                eidArrivalTimes[srcEid] = ns3::Simulator::Now();
            }

            packet->AddHeader(idp);
            auto res = bufferQueues[srcEid].insert(packet);
            if (res.second)
            {
                globalBufferQueue1 += Psize;
                maxGlobalBufferQueue1 = std::max(maxGlobalBufferQueue1, globalBufferQueue1);
            }
        }
        else
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
        }
    }

    else if (reorderAddr.GetIpv6() == Ipv6Address("2001:4::200:ff:fe00:7"))
    {
        Ptr<Packet> packet = sock->Recv();
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);
        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();
        uint8_t type = idp.GetType();
        static std::map<std::string, ns3::Time> eidArrivalTimes;

        if (expectedSeqNums.find(srcEid) == expectedSeqNums.end())
        {
            expectedSeqNums[srcEid] = 1;
        }

        if (seq == expectedSeqNums[srcEid])
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:4::200:ff:fe00:8"));
            expectedSeqNums[srcEid]++;
            while ((bufferQueues.find(srcEid) != bufferQueues.end()) && (!bufferQueues[srcEid].empty()))
            {
                auto pktIterator = bufferQueues[srcEid].begin();
                Ptr<Packet> bufferedPacket = *pktIterator;
                IDPHeader bufferedIdp;
                bufferedPacket->RemoveHeader(bufferedIdp);
                if (bufferedIdp.GetSeqNum() == expectedSeqNums[srcEid])
                {
                    uint32_t pktSize;
                    pktSize = bufferedIdp.GetPacketSize();
                    bufferedPacket->AddHeader(bufferedIdp);
                    sock->SendTo(bufferedPacket, 0, Inet6SocketAddress("2001:4::200:ff:fe00:8"));
                    expectedSeqNums[srcEid]++;
                    globalBufferQueue2 -= pktSize;
                    bufferQueues[srcEid].erase(pktIterator);

                    if (bufferQueues[srcEid].empty())
                    {
                        bufferQueues.erase(srcEid);
                    }
                }
                else
                {
                    bufferedPacket->AddHeader(bufferedIdp);
                    break;
                }
            }
        }

        else if (seq > expectedSeqNums[srcEid])
        {
            uint32_t Psize = idp.GetPacketSize();
            if (bufferQueues.find(srcEid) == bufferQueues.end())
            {
                bufferQueues[srcEid] = std::set<Ptr<Packet>, PacketCompare>();
                eidArrivalTimes[srcEid] = ns3::Simulator::Now();
            }

            packet->AddHeader(idp);
            auto res = bufferQueues[srcEid].insert(packet);
            if (res.second)
            {
                globalBufferQueue2 += Psize;
                maxGlobalBufferQueue2 = std::max(maxGlobalBufferQueue2, globalBufferQueue2);
            }
        }
        else
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:4::200:ff:fe00:8"));
        }
    }
}

std::map<std::string, int> Pk;
std::map<std::string, bool> Qk;
std::map<std::string, bool> Mk;
std::map<std::string, Time> Rk;

void recvwithMYRORDER(Ptr<Socket> sock)
{
    double minPkt = 64.0;
    double maxPkt = 1024.0;
    static Time last_time = Seconds(0);
    double p = 0.0;
    static double pktsize_thresh = 0.0;
    static double Lnorm = 0.0;
    static double gamma = 0.8;
    static double lambda = 0.0;
    static double curMinScore = 0.0;
    static uint32_t N_pl = 0;
    uint32_t N_rt = 0;
    static uint32_t standbyBuffer = 0.1 * maxBufferSize;
    static uint32_t globalBufferQueue;
    static std::map<std::string, std::set<Ptr<Packet>, PacketCompare>> bufferQueues;
    static std::map<std::string, uint32_t> expectedSeqNums;

    Address boundAddress;
    sock->GetSockName(boundAddress);
    Inet6SocketAddress reorderAddr = Inet6SocketAddress::ConvertFrom(boundAddress);
    if (reorderAddr.GetIpv6() == Ipv6Address("2001:3::200:ff:fe00:5"))
    {
        Ptr<Packet> packet = sock->Recv();
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);
        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();
        uint8_t type = idp.GetType();
        uint32_t N;

        if (expectedSeqNums.find(srcEid) == expectedSeqNums.end())
        {
            expectedSeqNums[srcEid] = 1;
        }

        if (seq == expectedSeqNums[srcEid])
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
            expectedSeqNums[srcEid]++;
            while ((bufferQueues.find(srcEid) != bufferQueues.end()) && (!bufferQueues[srcEid].empty()))
            {
                auto pktIterator = bufferQueues[srcEid].begin();
                Ptr<Packet> bufferedPacket = *pktIterator;
                IDPHeader bufferedIdp;
                bufferedPacket->RemoveHeader(bufferedIdp);
                if (bufferedIdp.GetSeqNum() == expectedSeqNums[srcEid])
                {
                    uint32_t pktSize;
                    pktSize = bufferedIdp.GetPacketSize();
                    bufferedPacket->AddHeader(bufferedIdp);

                    sock->SendTo(bufferedPacket, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
                    expectedSeqNums[srcEid]++;
                    globalBufferQueue -= pktSize;
                    bufferQueues[srcEid].erase(pktIterator);

                    if (bufferQueues[srcEid].empty())
                    {
                        Pk.erase(srcEid);
                        Qk.erase(srcEid);
                        Mk.erase(srcEid);
                        Rk.erase(srcEid);
                        for (auto it = scoreSet.begin(); it != scoreSet.end(); ++it)
                        {
                            if (it->srcEid == srcEid)
                            {
                                scoreSet.erase(it);
                                break;
                            }
                        }
                        bufferQueues.erase(srcEid);
                    }
                    else
                    {
                        Lnorm = 0.1 + ((static_cast<double>(pktSize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                        calScore(seq, srcEid, bufferQueues, expectedSeqNums, Lnorm, 1);
                    }
                }
                else
                {
                    bufferedPacket->AddHeader(bufferedIdp);
                    break;
                }
            }
        }

        else if (seq > expectedSeqNums[srcEid])
        {
            uint32_t Fsizecur;
            uint32_t Psize = idp.GetPacketSize();
            if (bufferQueues.find(srcEid) == bufferQueues.end())
            {
                bufferQueues[srcEid] = std::set<Ptr<Packet>, PacketCompare>();
                uint32_t Psize = idp.GetPacketSize();
                Pk[srcEid] = false;
                Qk[srcEid] = false;
                Mk[srcEid] = false;
                Rk[srcEid] = Time::Max();
                if (Psize <= pktsize_thresh)
                {
                    Qk[srcEid] = true;
                }
            }
            if (maxBufferSize <= globalBufferQueue)
            {
                cout << "full." << endl;
                cout << "Eid : " << srcEid << " 's " << idp.GetSeqNum() << "packet." << endl;
                packet->AddHeader(idp);
                int sendbyte = sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
                cout << "send: " << sendbyte << endl;
            }
            else
            {
                if (maxBufferSize - globalBufferQueue > standbyBuffer)
                {

                    packet->AddHeader(idp);
                    auto res = bufferQueues[srcEid].insert(packet);

                    if (res.second)
                    {
                        globalBufferQueue += Psize;
                        Lnorm = 0.1 + ((static_cast<double>(Psize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                        calScore(seq, srcEid, bufferQueues, expectedSeqNums, Lnorm, 0);
                    }
                }
                else if (maxBufferSize - globalBufferQueue <= standbyBuffer)
                {
                    packet->AddHeader(idp);
                    auto res = bufferQueues[srcEid].insert(packet);
                    if (res.second)
                    {
                        globalBufferQueue += Psize;
                        Lnorm = 0.1 + ((static_cast<double>(Psize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                        calScore(seq, srcEid, bufferQueues, expectedSeqNums, Lnorm, 0);
                    }
                    if (N_pl != 0)
                    {
                        string maxEid = getMaxScoreSrcEid();
                        uint32_t frdPsize;
                        if (!bufferQueues[maxEid].empty())
                        {
                            uint32_t totalSentSize = 0;
                            while (!bufferQueues[maxEid].empty() && totalSentSize < Psize)
                            {
                                auto forwardIterator = bufferQueues[maxEid].begin();
                                Ptr<Packet> forwardPacket = *forwardIterator;
                                IDPHeader forwardIdp;
                                forwardPacket->RemoveHeader(forwardIdp);
                                frdPsize = forwardIdp.GetPacketSize();
                                forwardPacket->AddHeader(forwardIdp);
                                if (totalSentSize > Psize)
                                {
                                    break;
                                }
                                sock->SendTo(forwardPacket, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
                                globalBufferQueue -= frdPsize;
                                totalSentSize += frdPsize;
                                bufferQueues[maxEid].erase(forwardIterator);

                                if (bufferQueues[maxEid].empty())
                                {
                                    Pk.erase(maxEid);
                                    Qk.erase(maxEid);
                                    Mk.erase(maxEid);
                                    Rk.erase(maxEid);
                                    for (auto it = scoreSet.begin(); it != scoreSet.end(); ++it)
                                    {
                                        if (it->srcEid == maxEid)
                                        {
                                            scoreSet.erase(it);
                                            break;
                                        }
                                    }
                                    bufferQueues.erase(maxEid);
                                    break;
                                }
                            }
                            if (!bufferQueues[maxEid].empty())
                            {
                                Lnorm = 0.1 + ((static_cast<double>(frdPsize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                                calScore(seq, maxEid, bufferQueues, expectedSeqNums, Lnorm, 1);
                            }
                        }
                    }

                    else
                    {
                        string maxEid = getMaxScoreSrcEid();
                        curMinScore = getScoreBySrcEid(getMinScoreSrcEid());
                        p = CalculateProbabilityForEid(srcEid, scoreSet);
                        std::random_device rd;
                        std::mt19937 gen(rd());
                        std::uniform_real_distribution<> distrib(0, 1);
                        if (distrib(gen) < p)
                        {
                            string maxEid = getMaxScoreSrcEid();
                            if (!bufferQueues[maxEid].empty())
                            {
                                uint32_t totalSentSize = 0;
                                uint32_t frdPsize;
                                while (!bufferQueues[maxEid].empty() && totalSentSize < Psize)
                                {
                                    auto forwardIterator = bufferQueues[maxEid].begin();
                                    Ptr<Packet> forwardPacket = *forwardIterator;
                                    IDPHeader forwardIdp;
                                    forwardPacket->RemoveHeader(forwardIdp);
                                    frdPsize = forwardIdp.GetPacketSize();
                                    forwardPacket->AddHeader(forwardIdp);
                                    if (totalSentSize > Psize)
                                    {
                                        break;
                                    }
                                    sock->SendTo(forwardPacket, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
                                    globalBufferQueue -= frdPsize;
                                    totalSentSize += frdPsize;
                                    bufferQueues[maxEid].erase(forwardIterator);
                                }
                                if (bufferQueues[maxEid].empty())
                                {
                                    Pk.erase(maxEid);
                                    Qk.erase(maxEid);
                                    Mk.erase(maxEid);
                                    Rk.erase(maxEid);
                                    for (auto it = scoreSet.begin(); it != scoreSet.end(); ++it)
                                    {
                                        if (it->srcEid == maxEid)
                                        {
                                            scoreSet.erase(it);
                                            break;
                                        }
                                    }
                                    bufferQueues.erase(maxEid);
                                }
                                else
                                {
                                    Lnorm = 0.1 + ((static_cast<double>(frdPsize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                                    calScore(seq, maxEid, bufferQueues, expectedSeqNums, Lnorm, 1);
                                }
                            }
                        }
                    }
                }
            }

            Time now = Simulator::Now();
            Fsizecur = bufferQueues[srcEid].size();
            Pk[srcEid] = IsPersistLargeFlow(srcEid, Fsizecur, Psize, standbyBuffer, Mk, Rk, now);
            if (!bufferQueues[srcEid].empty() && now > last_time + MilliSeconds(10))
            {
                pktsize_thresh = gamma * pktsize_thresh + (1 - gamma) * Psize;
                Fsizecur = bufferQueues[srcEid].size();
                N_pl = 0;
                for (const auto &entry : Pk)
                {
                    if (entry.second)
                    {
                        N_pl++;
                    }
                }
                for (const auto &entry : Qk)
                {
                    if (entry.second)
                    {
                        N_rt++;
                    }
                }
                N = bufferQueues.size();
                lambda = alpha * (static_cast<double>(N) / (1 + N_pl)) * sqrt(static_cast<double>(N_rt) / (N - N_rt + 1));
                standbyBuffer = gamma * standbyBuffer + (1 - gamma) * (maxBufferSize / (1 + lambda));
                last_time = now;
            }
        }
        else
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:3::200:ff:fe00:6"));
        }
    }

    else if (reorderAddr.GetIpv6() == Ipv6Address("2001:4::200:ff:fe00:7"))
    {
        Ptr<Packet> packet = sock->Recv();
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);
        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();
        uint8_t type = idp.GetType();
        uint32_t N;

        if (expectedSeqNums.find(srcEid) == expectedSeqNums.end())
        {
            expectedSeqNums[srcEid] = 1;
        }

        if (seq == expectedSeqNums[srcEid])
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:4::200:ff:fe00:8"));
            expectedSeqNums[srcEid]++;
            while ((bufferQueues.find(srcEid) != bufferQueues.end()) && (!bufferQueues[srcEid].empty()))
            {
                auto pktIterator = bufferQueues[srcEid].begin();
                Ptr<Packet> bufferedPacket = *pktIterator;
                IDPHeader bufferedIdp;
                bufferedPacket->RemoveHeader(bufferedIdp);
                if (bufferedIdp.GetSeqNum() == expectedSeqNums[srcEid])
                {
                    uint32_t pktSize;
                    pktSize = bufferedIdp.GetPacketSize();
                    bufferedPacket->AddHeader(bufferedIdp);
                    sock->SendTo(bufferedPacket, 0, Inet6SocketAddress("2001:4::200:ff:fe00:8"));

                    while (forwardedSeqNumsTable[srcEid].find(expectedSeqNums[srcEid] + 1) != forwardedSeqNumsTable[srcEid].end())
                    {
                        expectedSeqNums[srcEid]++;
                    }
                    expectedSeqNums[srcEid]++;
                    globalBufferQueue -= pktSize;
                    bufferQueues[srcEid].erase(pktIterator);

                    if (bufferQueues[srcEid].empty())
                    {
                        Pk.erase(srcEid);
                        Qk.erase(srcEid);
                        Mk.erase(srcEid);
                        Rk.erase(srcEid);
                        for (auto it = scoreSet.begin(); it != scoreSet.end(); ++it)
                        {
                            if (it->srcEid == srcEid)
                            {
                                scoreSet.erase(it);
                                break;
                            }
                        }
                        bufferQueues.erase(srcEid);
                    }
                    else
                    {
                        Lnorm = 0.1 + ((static_cast<double>(pktSize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                        calScore(seq, srcEid, bufferQueues, expectedSeqNums, Lnorm, 1);
                    }
                }
                else
                {
                    bufferedPacket->AddHeader(bufferedIdp);
                    break;
                }
            }
        }

        else if (seq > expectedSeqNums[srcEid])
        {
            uint32_t Fsizecur;
            uint32_t Psize = idp.GetPacketSize();
            if (bufferQueues.find(srcEid) == bufferQueues.end())
            {

                bufferQueues[srcEid] = std::set<Ptr<Packet>, PacketCompare>();
                uint32_t Psize = idp.GetPacketSize();
                Pk[srcEid] = false;
                Qk[srcEid] = false;
                Mk[srcEid] = false;
                Rk[srcEid] = Time::Max();
                if (Psize <= pktsize_thresh)
                {
                    Qk[srcEid] = true;
                }
            }
            if (maxBufferSize <= globalBufferQueue)
            {
                cout << "full." << endl;
                cout << "Eid : " << srcEid << " 's " << idp.GetSeqNum() << "packet." << endl;
                packet->AddHeader(idp);
                int sendbyte = sock->SendTo(packet, 0, Inet6SocketAddress("2001:4::200:ff:fe00:8"));
                cout << "send: " << sendbyte << endl;
            }
            else
            {
                if (maxBufferSize - globalBufferQueue > standbyBuffer)
                {

                    packet->AddHeader(idp);
                    auto res = bufferQueues[srcEid].insert(packet);

                    if (res.second)
                    {
                        globalBufferQueue += Psize;
                        Lnorm = 0.1 + ((static_cast<double>(Psize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                        calScore(seq, srcEid, bufferQueues, expectedSeqNums, Lnorm, 0);
                    }
                }
                else if (maxBufferSize - globalBufferQueue <= standbyBuffer)
                {
                    packet->AddHeader(idp);
                    auto res = bufferQueues[srcEid].insert(packet);
                    if (res.second)
                    {
                        globalBufferQueue += Psize;
                        Lnorm = 0.1 + ((static_cast<double>(Psize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                        calScore(seq, srcEid, bufferQueues, expectedSeqNums, Lnorm, 0);
                    }

                    if (N_pl != 0)
                    {
                        string maxEid = getMaxScoreSrcEid();
                        uint32_t frdPsize;
                        if (!bufferQueues[maxEid].empty())
                        {
                            uint32_t totalSentSize = 0;
                            while (!bufferQueues[maxEid].empty() && totalSentSize < Psize)
                            {
                                auto forwardIterator = bufferQueues[maxEid].begin();
                                Ptr<Packet> forwardPacket = *forwardIterator;
                                IDPHeader forwardIdp;
                                forwardPacket->RemoveHeader(forwardIdp);
                                frdPsize = forwardIdp.GetPacketSize();
                                forwardPacket->AddHeader(forwardIdp);
                                if (totalSentSize > Psize)
                                {
                                    break;
                                }
                                sock->SendTo(forwardPacket, 0, Inet6SocketAddress("2001:4::200:ff:fe00:8"));

                                globalBufferQueue -= frdPsize;
                                totalSentSize += frdPsize;
                                bufferQueues[maxEid].erase(forwardIterator);

                                if (bufferQueues[maxEid].empty())
                                {
                                    Pk.erase(maxEid);
                                    Qk.erase(maxEid);
                                    Mk.erase(maxEid);
                                    Rk.erase(maxEid);
                                    for (auto it = scoreSet.begin(); it != scoreSet.end(); ++it)
                                    {
                                        if (it->srcEid == maxEid)
                                        {
                                            scoreSet.erase(it);
                                            break;
                                        }
                                    }
                                    bufferQueues.erase(maxEid);
                                    break;
                                }
                            }

                            if (!bufferQueues[maxEid].empty())
                            {
                                Lnorm = 0.1 + ((static_cast<double>(frdPsize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                                calScore(seq, maxEid, bufferQueues, expectedSeqNums, Lnorm, 1);
                            }
                        }
                    }

                    else
                    {
                        string maxEid = getMaxScoreSrcEid();
                        curMinScore = getScoreBySrcEid(getMinScoreSrcEid());
                        p = CalculateProbabilityForEid(srcEid, scoreSet);
                        std::random_device rd;
                        std::mt19937 gen(rd());
                        std::uniform_real_distribution<> distrib(0, 1);
                        if (distrib(gen) < p)
                        {
                            string maxEid = getMaxScoreSrcEid();
                            if (!bufferQueues[maxEid].empty())
                            {
                                uint32_t totalSentSize = 0;
                                uint32_t frdPsize;
                                while (!bufferQueues[maxEid].empty() && totalSentSize < Psize)
                                {
                                    auto forwardIterator = bufferQueues[maxEid].begin();
                                    Ptr<Packet> forwardPacket = *forwardIterator;
                                    IDPHeader forwardIdp;
                                    forwardPacket->RemoveHeader(forwardIdp);
                                    frdPsize = forwardIdp.GetPacketSize();
                                    forwardPacket->AddHeader(forwardIdp);
                                    if (totalSentSize > Psize)
                                    {
                                        break;
                                    }

                                    sock->SendTo(forwardPacket, 0, Inet6SocketAddress("2001:4::200:ff:fe00:8"));
                                    globalBufferQueue -= frdPsize;
                                    totalSentSize += frdPsize;

                                    bufferQueues[maxEid].erase(forwardIterator);
                                }
                                if (bufferQueues[maxEid].empty())
                                {

                                    Pk.erase(maxEid);
                                    Qk.erase(maxEid);
                                    Mk.erase(maxEid);
                                    Rk.erase(maxEid);
                                    for (auto it = scoreSet.begin(); it != scoreSet.end(); ++it)
                                    {
                                        if (it->srcEid == maxEid)
                                        {
                                            scoreSet.erase(it);
                                            break;
                                        }
                                    }
                                    bufferQueues.erase(maxEid);
                                }
                                else
                                {
                                    Lnorm = 0.1 + ((static_cast<double>(frdPsize) - minPkt) / (maxPkt - minPkt)) * 0.9;
                                    calScore(seq, maxEid, bufferQueues, expectedSeqNums, Lnorm, 1);
                                }
                            }
                        }
                    }
                }

                else
                {
                    cout << "maxBufferSize: " << maxBufferSize << endl;
                    cout << "globalBufferQueue" << globalBufferQueue << endl;
                    cout << "Standbybuffer" << standbyBuffer << endl;
                }
            }
            Time now = Simulator::Now();
            Fsizecur = bufferQueues[srcEid].size();
            Pk[srcEid] = IsPersistLargeFlow(srcEid, Fsizecur, Psize, standbyBuffer, Mk, Rk, now);
            if (!bufferQueues[srcEid].empty() && now > last_time + MilliSeconds(10))
            {
                pktsize_thresh = gamma * pktsize_thresh + (1 - gamma) * Psize;
                Fsizecur = bufferQueues[srcEid].size();
                N_pl = 0;
                for (const auto &entry : Pk)
                {
                    if (entry.second)
                    {
                        N_pl++;
                    }
                }
                for (const auto &entry : Qk)
                {
                    if (entry.second)
                    {
                        N_rt++;
                    }
                }
                N = bufferQueues.size();
                lambda = alpha * (static_cast<double>(N) / (1 + N_pl)) * sqrt(static_cast<double>(N_rt) / (N - N_rt + 1));
                last_time = now;
            }
        }
        else
        {
            packet->AddHeader(idp);
            sock->SendTo(packet, 0, Inet6SocketAddress("2001:4::200:ff:fe00:8"));
        }
    }
}

void recv2(Ptr<Socket> sock)
{

    static double receivedBytes = 0;
    Time firstReceiveTime = Simulator::Now();
    static Time lastReceiveTime = Seconds(0);
    if (firstReceiveTime > lastReceiveTime + Seconds(0.01))
    {

        double duration = (firstReceiveTime - lastReceiveTime).GetSeconds();
        double throughput = (receivedBytes * 8) / (duration * 1e6);
        outfile3 << firstReceiveTime.GetSeconds() << "," << throughput << "\n";
        receivedBytes = 0;
        lastReceiveTime = firstReceiveTime;
    }
    static std::map<std::string, uint32_t> expectedSeqNums;
    static std::map<std::string, std::set<Ptr<Packet>, PacketCompare>> RbufferQueues;

    Address boundAddress;
    sock->GetSockName(boundAddress);
    Inet6SocketAddress reorderAddr = Inet6SocketAddress::ConvertFrom(boundAddress);
    if (reorderAddr.GetIpv6() == Ipv6Address("2001:3::200:ff:fe00:6"))
    {
        Ptr<Packet> packet = sock->Recv();
        if (packet == nullptr)
        {
            cout << "Received a null packet pointer." << endl;
            return;
        }
        // totalNum++;
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);

        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();

        if (expectedSeqNums.find(srcEid) == expectedSeqNums.end())
        {
            expectedSeqNums[srcEid] = 1;
        }

        if (seq == expectedSeqNums[srcEid])
        {
            receivedBytes += idp.GetPacketSize();
            expectedSeqNums[srcEid]++;
            while (!RbufferQueues[srcEid].empty())
            {
                auto pktIterator = RbufferQueues[srcEid].begin();
                Ptr<Packet> bufferedPacket = *pktIterator;
                IDPHeader bufferedIdp;
                bufferedPacket->RemoveHeader(bufferedIdp);
                if (bufferedIdp.GetSeqNum() == expectedSeqNums[srcEid])
                {
                    // ②
                    receivedBytes += bufferedIdp.GetPacketSize();
                    expectedSeqNums[srcEid]++;
                    Ptr<Packet> msg_packet = Create<Packet>(52);
                    Ipv6Header ipv6Header;
                    ipv6Header.SetSourceAddress(Ipv6Address("2001:3::200:ff:fe00:6"));
                    ipv6Header.SetDestinationAddress(Ipv6Address("2001:3::200:ff:fe00:5"));
                    ipv6Header.SetPayloadLength(packet->GetSize());
                    msg_packet->AddHeader(ipv6Header);
                    IDPHeader idpHeader;
                    idpHeader.SetDestinationEid(recveid_0);
                    idpHeader.SetType('M');
                    idpHeader.SetSourceEid(srcEid);
                    idpHeader.SetSeqNum(bufferedIdp.GetSeqNum());
                    msg_packet->AddHeader(idpHeader);
                    Ptr<Node> n = sock->GetNode();
                    Ptr<Socket> msgSock = Socket::CreateSocket(n, Ipv6RawSocketFactory::GetTypeId());
                    msgSock->Bind();
                    msgSock->Connect(Inet6SocketAddress("2001:3::200:ff:fe00:5"));
                    msgSock->SetAttribute("Protocol", UintegerValue(153));
                    msgSock->Send(msg_packet);
                    msgSock->Close();
                    RbufferQueues[srcEid].erase(pktIterator);
                }
                else
                {
                    bufferedPacket->AddHeader(bufferedIdp);
                    break;
                }
            }
            if (RbufferQueues[srcEid].size() > 0)
            {
                auto minSeqPkt = *RbufferQueues[srcEid].begin();
                IDPHeader minSeqPktIdp;
                minSeqPkt->RemoveHeader(minSeqPktIdp);
                uint32_t minSeqNum = minSeqPktIdp.GetSeqNum();

                if (minSeqNum > expectedSeqNums[srcEid] + 750)
                {
                    minSeqPkt->AddHeader(minSeqPktIdp);
                    for (auto &resendPacket : RbufferQueues[srcEid])
                    {
                        IDPHeader resendIdp;
                        resendPacket->RemoveHeader(resendIdp);
                        resendIdp.SetDestinationEid(sendeid_0);
                        resendIdp.SetType('R');
                        uint32_t Rseq = resendIdp.GetSeqNum();
                        uint32_t headersSize = resendIdp.GetSerializedSize();
                        uint32_t packetSize = resendPacket->GetSize();
                        if (packetSize > headersSize)
                        {
                            resendPacket->RemoveAtEnd(packetSize - headersSize);
                        }
                        cout << "resend request for Eid " << srcEid << " No. " << Rseq << endl;
                        resendPacket->AddHeader(resendIdp);

                        std::random_device rd;
                        std::mt19937 gen(rd());
                        std::uniform_int_distribution<> distrib(0, 9);
                        int random_number = distrib(gen);
                        Ptr<Node> n = sock->GetNode();
                        Ptr<Socket> retransSock = Socket::CreateSocket(n, Ipv6RawSocketFactory::GetTypeId());
                        retransSock->Bind();
                        retransSock->Connect(RelayAddresses[random_number]);
                        retransSock->SetAttribute("Protocol", UintegerValue(153));
                        retransSock->Send(resendPacket);
                        retransSock->Close();
                    }
                    RbufferQueues[srcEid].clear();
                    return;
                }
                else
                {
                    minSeqPkt->AddHeader(minSeqPktIdp);
                    return;
                }
            }
        }
        else if (seq > expectedSeqNums[srcEid])
        {
            packet->AddHeader(idp);
            RbufferQueues[srcEid].insert(packet);

            if (RbufferQueues[srcEid].size() > 0)
            {
                auto minSeqPkt = *RbufferQueues[srcEid].begin();
                IDPHeader minSeqPktIdp;
                minSeqPkt->RemoveHeader(minSeqPktIdp);
                uint32_t minSeqNum = minSeqPktIdp.GetSeqNum();

                if (minSeqNum > expectedSeqNums[srcEid] + 5)
                {
                    minSeqPkt->AddHeader(minSeqPktIdp);
                    for (auto &resendPacket : RbufferQueues[srcEid])
                    {
                        IDPHeader resendIdp;
                        resendPacket->RemoveHeader(resendIdp);
                        resendIdp.SetDestinationEid(sendeid_0);
                        uint32_t Rseq = resendIdp.GetSeqNum();
                        resendIdp.SetType('R');
                        cout << "resend request for Eid " << srcEid << " No. " << Rseq << endl;
                        uint32_t headersSize = resendIdp.GetSerializedSize();
                        uint32_t packetSize = resendPacket->GetSize();
                        if (packetSize > headersSize)
                        {
                            resendPacket->RemoveAtEnd(packetSize - headersSize);
                        }
                        resendPacket->AddHeader(resendIdp);
                        std::random_device rd;
                        std::mt19937 gen(rd());
                        std::uniform_int_distribution<> distrib(0, 9);
                        int random_number = distrib(gen);
                        Ptr<Node> n = sock->GetNode();
                        Ptr<Socket> retransSock = Socket::CreateSocket(n, Ipv6RawSocketFactory::GetTypeId());
                        retransSock->Bind();
                        retransSock->Connect(RelayAddresses[random_number]);
                        retransSock->SetAttribute("Protocol", UintegerValue(153));
                        retransSock->Send(resendPacket);
                        retransSock->Close();
                    }
                    RbufferQueues[srcEid].clear();
                    return;
                }
                else
                {
                    minSeqPkt->AddHeader(minSeqPktIdp);
                    return;
                }
            }
        }
        else
        {
            return;
        }
    }
    else if (reorderAddr.GetIpv6() == Ipv6Address("2001:4::200:ff:fe00:8"))
    {
        Ptr<Packet> packet = sock->Recv();
        Ipv6Header ipHdr;
        packet->RemoveHeader(ipHdr);
        IDPHeader idp;
        packet->RemoveHeader(idp);

        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        uint32_t seq = idp.GetSeqNum();

        if (expectedSeqNums.find(srcEid) == expectedSeqNums.end())
        {
            expectedSeqNums[srcEid] = 1;
        }

        if (seq == expectedSeqNums[srcEid])
        {
            receivedBytes += idp.GetPacketSize();
            expectedSeqNums[srcEid]++;
            while (!RbufferQueues[srcEid].empty())
            {
                auto pktIterator = RbufferQueues[srcEid].begin();
                Ptr<Packet> bufferedPacket = *pktIterator;
                IDPHeader bufferedIdp;
                bufferedPacket->RemoveHeader(bufferedIdp);
                if (bufferedIdp.GetSeqNum() == expectedSeqNums[srcEid])
                {
                    // ②
                    receivedBytes += bufferedIdp.GetPacketSize();
                    expectedSeqNums[srcEid]++;
                    Ptr<Packet> msg_packet = Create<Packet>(52);
                    Ipv6Header ipv6Header;
                    ipv6Header.SetSourceAddress(Ipv6Address("2001:4::200:ff:fe00:8"));
                    ipv6Header.SetDestinationAddress(Ipv6Address("2001:4::200:ff:fe00:7"));
                    ipv6Header.SetPayloadLength(packet->GetSize());
                    msg_packet->AddHeader(ipv6Header);
                    IDPHeader idpHeader;
                    idpHeader.SetDestinationEid(recveid_0);
                    idpHeader.SetType('M');
                    idpHeader.SetSourceEid(srcEid);
                    idpHeader.SetSeqNum(bufferedIdp.GetSeqNum());
                    msg_packet->AddHeader(idpHeader);
                    Ptr<Node> n = sock->GetNode();
                    Ptr<Socket> msgSock = Socket::CreateSocket(n, Ipv6RawSocketFactory::GetTypeId());
                    msgSock->Bind();
                    msgSock->Connect(Inet6SocketAddress("2001:4::200:ff:fe00:7"));
                    msgSock->SetAttribute("Protocol", UintegerValue(153));
                    msgSock->Send(msg_packet);

                    msgSock->Close();
                    RbufferQueues[srcEid].erase(pktIterator);
                }
                else
                {
                    bufferedPacket->AddHeader(bufferedIdp);
                    break;
                }
            }
            if (RbufferQueues[srcEid].size() > 0)
            {
                auto minSeqPkt = *RbufferQueues[srcEid].begin();
                IDPHeader minSeqPktIdp;
                minSeqPkt->RemoveHeader(minSeqPktIdp);
                uint32_t minSeqNum = minSeqPktIdp.GetSeqNum();

                if (minSeqNum > expectedSeqNums[srcEid] + 5)
                {
                    minSeqPkt->AddHeader(minSeqPktIdp);
                    for (auto &resendPacket : RbufferQueues[srcEid])
                    {
                        IDPHeader resendIdp;
                        resendPacket->RemoveHeader(resendIdp);
                        resendIdp.SetDestinationEid(sendeid_1);
                        resendIdp.SetType('R');
                        uint32_t Rseq = resendIdp.GetSeqNum();
                        uint32_t headersSize = resendIdp.GetSerializedSize();
                        uint32_t packetSize = resendPacket->GetSize();
                        if (packetSize > headersSize)
                        {
                            resendPacket->RemoveAtEnd(packetSize - headersSize);
                        }
                        cout << "resend request for Eid " << srcEid << " No. " << Rseq << endl;
                        resendPacket->AddHeader(resendIdp);

                        std::random_device rd;
                        std::mt19937 gen(rd());
                        std::uniform_int_distribution<> distrib(0, 9);
                        int random_number = distrib(gen);
                        Ptr<Node> n = sock->GetNode();
                        Ptr<Socket> retransSock = Socket::CreateSocket(n, Ipv6RawSocketFactory::GetTypeId());
                        retransSock->Bind();
                        retransSock->Connect(RelayAddresses[random_number]);
                        retransSock->SetAttribute("Protocol", UintegerValue(153));
                        retransSock->Send(resendPacket);
                        retransSock->Close();
                    }
                    RbufferQueues[srcEid].clear();
                    return;
                }
                else
                {
                    minSeqPkt->AddHeader(minSeqPktIdp);
                    return;
                }
            }
        }
        else if (seq > expectedSeqNums[srcEid])
        {
            packet->AddHeader(idp);
            RbufferQueues[srcEid].insert(packet);

            if (RbufferQueues[srcEid].size() > 0)
            {
                auto minSeqPkt = *RbufferQueues[srcEid].begin();
                IDPHeader minSeqPktIdp;
                minSeqPkt->RemoveHeader(minSeqPktIdp);
                uint32_t minSeqNum = minSeqPktIdp.GetSeqNum();

                if (minSeqNum > expectedSeqNums[srcEid] + 750)
                {
                    minSeqPkt->AddHeader(minSeqPktIdp);
                    for (auto &resendPacket : RbufferQueues[srcEid])
                    {
                        IDPHeader resendIdp;
                        resendPacket->RemoveHeader(resendIdp);
                        resendIdp.SetDestinationEid(sendeid_1);
                        uint32_t Rseq = resendIdp.GetSeqNum();
                        resendIdp.SetType('R');
                        uint32_t headersSize = resendIdp.GetSerializedSize();
                        uint32_t packetSize = resendPacket->GetSize();
                        if (packetSize > headersSize)
                        {
                            resendPacket->RemoveAtEnd(packetSize - headersSize);
                            resendPacket->AddHeader(resendIdp);
                            std::random_device rd;
                            std::mt19937 gen(rd());
                            std::uniform_int_distribution<> distrib(0, 9);
                            int random_number = distrib(gen);
                            Ptr<Node> n = sock->GetNode();
                            Ptr<Socket> retransSock = Socket::CreateSocket(n, Ipv6RawSocketFactory::GetTypeId());
                            retransSock->Bind();
                            retransSock->Connect(RelayAddresses[random_number]);
                            retransSock->SetAttribute("Protocol", UintegerValue(153));
                            retransSock->Send(resendPacket);
                            retransSock->Close();
                        }
                        RbufferQueues[srcEid].clear();
                        return;
                    }
                    else
                    {
                        minSeqPkt->AddHeader(minSeqPktIdp);
                        return;
                    }
                }
            }
            else
            {
                return;
            }
        }
        else
        {
            return;
        }
    }
}
void resend(Ptr<Socket> sock)
{
    Address boundAddress;
    sock->GetSockName(boundAddress);
    Inet6SocketAddress recvAddr = Inet6SocketAddress::ConvertFrom(boundAddress);
    Ptr<Packet> packet = sock->Recv();
    Ipv6Header ipHdr;
    packet->RemoveHeader(ipHdr);
    IDPHeader idp;
    packet->RemoveHeader(idp);

    if (recvAddr.GetIpv6() == Ipv6Address("2001:1::200:ff:fe00:1") && idp.GetType() == 'R')
    {

        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        idp.SetDestinationEid(recveid_0);
        idp.SetType('D');
        packet->AddHeader(idp);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(0, 9);
        int random_number = distrib(gen);
        Ptr<Node> n = sock->GetNode();
        Ptr<Socket> retransSock = Socket::CreateSocket(n, Ipv6RawSocketFactory::GetTypeId());
        retransSock->Bind();
        retransSock->Connect(RelayAddresses[random_number]);
        retransSock->SetAttribute("Protocol", UintegerValue(153));
        uint32_t bytesSent = retransSock->Send(packet);

        retransSock->Close();
    }
    else if (recvAddr.GetIpv6() == Ipv6Address("2001:2::200:ff:fe00:3") && idp.GetType() == 'R')
    {
        std::string srcEid = idp.GetSourceEid().substr(0, 20);
        idp.SetDestinationEid(recveid_1);
        idp.SetType('D');
        packet->AddHeader(idp);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(0, 9);
        int random_number = distrib(gen);
        Ptr<Node> n = sock->GetNode();
        Ptr<Socket> retransSock = Socket::CreateSocket(n, Ipv6RawSocketFactory::GetTypeId());
        retransSock->Bind();
        retransSock->Connect(RelayAddresses[random_number]);
        retransSock->SetAttribute("Protocol", UintegerValue(153));
        uint32_t bytesSent = retransSock->Send(packet);
        retransSock->Close();
    }
    else
    {
        cout << "no match dest Addr!" << endl;
    }
}

NS_LOG_COMPONENT_DEFINE("GEANTDEMO");

NodeContainer hosts, edges, routers;

int main(int argc, char *argv[])
{
    LogComponentEnable("GEANTDEMO", LOG_LEVEL_ALL);

    NS_LOG_INFO("Create nodes.");
    std::string mode = "NOREORDER";
    bool pcap = false;
    outfile3 << "Time(s),Throughput\n";

    hosts.Create(4);
    edges.Create(4);
    routers.Create(28);

    std::vector<NodeContainer> nodeAdjancyList;
    // 0-3
    nodeAdjancyList.push_back(NodeContainer(hosts.Get(0), edges.Get(0))); // S1 - 1
    nodeAdjancyList.push_back(NodeContainer(hosts.Get(1), edges.Get(1))); // S2 - 4
    nodeAdjancyList.push_back(NodeContainer(edges.Get(2), hosts.Get(2))); // C1 - 25
    nodeAdjancyList.push_back(NodeContainer(edges.Get(3), hosts.Get(3))); // C2 - 29
    // 4-7
    nodeAdjancyList.push_back(NodeContainer(edges.Get(0), routers.Get(0))); // 1 - 2
    nodeAdjancyList.push_back(NodeContainer(edges.Get(0), routers.Get(1))); // 1 - 3
    nodeAdjancyList.push_back(NodeContainer(edges.Get(1), routers.Get(1))); // 4 - 3
    nodeAdjancyList.push_back(NodeContainer(edges.Get(1), routers.Get(2))); // 4 - 5
    // 8-12
    nodeAdjancyList.push_back(NodeContainer(routers.Get(0), routers.Get(1))); // 2 - 3
    nodeAdjancyList.push_back(NodeContainer(routers.Get(0), routers.Get(6))); // 2 - 9 √ 2
    nodeAdjancyList.push_back(NodeContainer(routers.Get(0), routers.Get(7))); // 2 - 10
    nodeAdjancyList.push_back(NodeContainer(routers.Get(1), routers.Get(2))); // 3 - 5
    nodeAdjancyList.push_back(NodeContainer(routers.Get(1), routers.Get(8))); // 3 - 11
    // 13-17
    nodeAdjancyList.push_back(NodeContainer(routers.Get(2), routers.Get(3)));  // 5 - 6
    nodeAdjancyList.push_back(NodeContainer(routers.Get(2), routers.Get(4)));  // 5 - 7
    nodeAdjancyList.push_back(NodeContainer(routers.Get(2), routers.Get(13))); // 5 - 16
    nodeAdjancyList.push_back(NodeContainer(routers.Get(3), routers.Get(10))); // 6 - 13 √ 6
    nodeAdjancyList.push_back(NodeContainer(routers.Get(3), routers.Get(14))); // 6 - 17
    // 18-22
    nodeAdjancyList.push_back(NodeContainer(routers.Get(4), routers.Get(6)));  // 7 - 9
    nodeAdjancyList.push_back(NodeContainer(routers.Get(4), routers.Get(10))); // 7 - 13
    nodeAdjancyList.push_back(NodeContainer(routers.Get(5), routers.Get(6)));  // 8 - 9
    nodeAdjancyList.push_back(NodeContainer(routers.Get(5), routers.Get(8)));  // 8 - 11
    nodeAdjancyList.push_back(NodeContainer(routers.Get(6), routers.Get(11))); // 9 - 14
    // 23-27
    nodeAdjancyList.push_back(NodeContainer(routers.Get(7), routers.Get(8)));  // 10 - 11 √ 10
    nodeAdjancyList.push_back(NodeContainer(routers.Get(7), routers.Get(15))); // 10 - 18
    nodeAdjancyList.push_back(NodeContainer(routers.Get(8), routers.Get(9)));  // 11 - 12 √ 11
    nodeAdjancyList.push_back(NodeContainer(routers.Get(8), routers.Get(11))); // 11 - 14
    nodeAdjancyList.push_back(NodeContainer(routers.Get(8), routers.Get(16))); // 11 - 19
    // 28-32
    nodeAdjancyList.push_back(NodeContainer(routers.Get(9), routers.Get(12)));  // 12 - 15
    nodeAdjancyList.push_back(NodeContainer(routers.Get(10), routers.Get(11))); // 13 - 14  √ 13
    nodeAdjancyList.push_back(NodeContainer(routers.Get(10), routers.Get(14))); // 13 - 17
    nodeAdjancyList.push_back(NodeContainer(routers.Get(11), routers.Get(12))); // 14 - 15  √ 14
    nodeAdjancyList.push_back(NodeContainer(routers.Get(11), routers.Get(13))); // 14 - 16
    // 33-37
    nodeAdjancyList.push_back(NodeContainer(routers.Get(11), routers.Get(22))); // 14 - 26
    nodeAdjancyList.push_back(NodeContainer(routers.Get(12), routers.Get(17))); // 15 - 20  √ 15
    nodeAdjancyList.push_back(NodeContainer(routers.Get(13), routers.Get(14))); // 16 - 17
    nodeAdjancyList.push_back(NodeContainer(routers.Get(13), routers.Get(26))); // 16 - 31  √ 16
    nodeAdjancyList.push_back(NodeContainer(routers.Get(15), routers.Get(16))); // 18 - 19
    // 38-42
    nodeAdjancyList.push_back(NodeContainer(routers.Get(15), routers.Get(18))); // 18 - 21
    nodeAdjancyList.push_back(NodeContainer(routers.Get(16), routers.Get(17))); // 19 - 20
    nodeAdjancyList.push_back(NodeContainer(routers.Get(16), routers.Get(21))); // 19 - 24  √ 19
    nodeAdjancyList.push_back(NodeContainer(routers.Get(17), routers.Get(22))); // 20 - 26
    nodeAdjancyList.push_back(NodeContainer(routers.Get(18), routers.Get(19))); // 21 - 22
    // 43-47
    nodeAdjancyList.push_back(NodeContainer(routers.Get(19), routers.Get(20))); // 22 - 23
    nodeAdjancyList.push_back(NodeContainer(routers.Get(19), routers.Get(21))); // 22 - 24
    nodeAdjancyList.push_back(NodeContainer(routers.Get(21), routers.Get(23))); // 24 - 27
    nodeAdjancyList.push_back(NodeContainer(routers.Get(22), routers.Get(23))); // 26 - 27 √ 26
    nodeAdjancyList.push_back(NodeContainer(routers.Get(22), routers.Get(26))); // 26 - 31
    // 48-52
    nodeAdjancyList.push_back(NodeContainer(routers.Get(23), routers.Get(24))); // 27 - 28
    nodeAdjancyList.push_back(NodeContainer(routers.Get(24), routers.Get(25))); // 28 - 30
    nodeAdjancyList.push_back(NodeContainer(routers.Get(24), routers.Get(26))); // 28 - 31
    nodeAdjancyList.push_back(NodeContainer(routers.Get(25), routers.Get(27))); // 30 - 32
    nodeAdjancyList.push_back(NodeContainer(routers.Get(26), routers.Get(27))); // 31 - 32
    // 53-58
    nodeAdjancyList.push_back(NodeContainer(routers.Get(20), edges.Get(2)));
    nodeAdjancyList.push_back(NodeContainer(routers.Get(21), edges.Get(2)));
    nodeAdjancyList.push_back(NodeContainer(routers.Get(23), edges.Get(2)));
    nodeAdjancyList.push_back(NodeContainer(routers.Get(24), edges.Get(3)));
    nodeAdjancyList.push_back(NodeContainer(edges.Get(2), edges.Get(3)));

    NS_LOG_INFO("Create Channel");
    std::vector<PointToPointHelper> p2p(nodeAdjancyList.size());

    for (size_t i = 0; i < 4; i++)
    {
        p2p[i].SetDeviceAttribute("DataRate", StringValue("200Mbps"));
        p2p[i].SetChannelAttribute("Delay", StringValue("1ms"));
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(25.0, 35.0);
    std::uniform_real_distribution<> disDelay(1.0, 2.0);

    for (size_t i = 4; i < p2p.size(); i++)
    {
        double randomBandwidth = dis(gen);
        std::string bandwidthStr = std::to_string(randomBandwidth) + "Mbps";
        p2p[i].SetDeviceAttribute("DataRate", StringValue(bandwidthStr));
        double randomDelay = disDelay(gen);
        std::string delayStr = std::to_string(randomDelay) + "ms";
        p2p[i].SetChannelAttribute("Delay", StringValue(delayStr));
    }

    std::vector<NetDeviceContainer> devices(nodeAdjancyList.size());
    for (uint32_t i = 0; i < nodeAdjancyList.size(); i++)
    {
        devices[i] = p2p[i].Install(nodeAdjancyList[i]);
    }

    RipNgHelper ripNgRouting;
    Ipv6ListRoutingHelper listRh;
    listRh.Add(ripNgRouting, 0);
    Ipv6StaticRoutingHelper staticRh;
    listRh.Add(staticRh, 5);

    InternetStackHelper internet;
    internet.SetRoutingHelper(listRh);
    internet.Install(hosts);
    internet.Install(edges);
    internet.Install(routers);

    NS_LOG_INFO("Assign Address.");
    std::vector<Ipv6InterfaceContainer> interfaces(nodeAdjancyList.size());
    Ipv6AddressHelper ipv6;
    for (uint32_t i = 0; i < nodeAdjancyList.size(); i++)
    {
        std::ostringstream subset;
        subset << "2001:" << i + 1 << "::";
        ipv6.SetBase(subset.str().c_str(), Ipv6Prefix(64));
        interfaces[i] = ipv6.Assign(devices[i]);
    }

    for (uint32_t i = 0; i < interfaces.size(); i++)
    {
        for (uint32_t j = 0; j < interfaces[i].GetN(); j++)
        {
            NS_LOG_INFO("Interface " << i << " Address " << j << ": " << interfaces[i].GetAddress(j, 1));
        }
    }

    eidRegister(sendeid_0, interfaces[0].GetAddress(0, 1)); // 2001:1::200:ff:fe00:1 20000000000000000000
    eidRegister(sendeid_1, interfaces[1].GetAddress(0, 1)); // 2001:2::200:ff:fe00:3 20000000000000000001

    eidRegister(mpsid, interfaces[9].GetAddress(0, 1));  // 2001:10::200:ff:fe00:13
    eidRegister(mpsid, interfaces[16].GetAddress(0, 1)); // 2001:17::200:ff:fe00:21
    eidRegister(mpsid, interfaces[23].GetAddress(0, 1)); // 2001:24::200:ff:fe00:2f
    eidRegister(mpsid, interfaces[25].GetAddress(0, 1)); // 2001:26::200:ff:fe00:33
    eidRegister(mpsid, interfaces[29].GetAddress(0, 1)); // 2001:30::200:ff:fe00:3b
    eidRegister(mpsid, interfaces[31].GetAddress(0, 1)); // 2001:32::200:ff:fe00:3f
    eidRegister(mpsid, interfaces[34].GetAddress(0, 1)); // 2001:35::200:ff:fe00:45
    eidRegister(mpsid, interfaces[36].GetAddress(0, 1)); // 2001:37::200:ff:fe00:49
    eidRegister(mpsid, interfaces[40].GetAddress(0, 1)); // 2001:41::200:ff:fe00:51
    eidRegister(mpsid, interfaces[46].GetAddress(0, 1)); // 2001:47::200:ff:fe00:5d

    eidRegister(recveid_0, interfaces[2].GetAddress(0, 1)); // 2001:3::200:ff:fe00:5 10000000000000000000
    eidRegister(recveid_1, interfaces[3].GetAddress(0, 1)); // 2001:4::200:ff:fe00:7 10000000000000000001

    Ptr<Socket> relay1Sock;
    Address relay1Addr(Inet6SocketAddress(interfaces[9].GetAddress(0, 1)));
    relay1Sock = Socket::CreateSocket(routers.Get(0), Ipv6RawSocketFactory::GetTypeId());
    relay1Sock->Bind(relay1Addr);
    relay1Sock->SetAttribute("Protocol", UintegerValue(153));
    relay1Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> relay2Sock;
    Address relay2Addr(Inet6SocketAddress(interfaces[16].GetAddress(0, 1)));
    relay2Sock = Socket::CreateSocket(routers.Get(3), Ipv6RawSocketFactory::GetTypeId());
    relay2Sock->Bind(relay2Addr);
    relay2Sock->SetAttribute("Protocol", UintegerValue(153));
    relay2Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> relay3Sock;
    Address relay3Addr(Inet6SocketAddress(interfaces[23].GetAddress(0, 1)));
    relay3Sock = Socket::CreateSocket(routers.Get(7), Ipv6RawSocketFactory::GetTypeId());
    relay3Sock->Bind(relay3Addr);
    relay3Sock->SetAttribute("Protocol", UintegerValue(153));
    relay3Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> relay4Sock;
    Address relay4Addr(Inet6SocketAddress(interfaces[25].GetAddress(0, 1)));
    relay4Sock = Socket::CreateSocket(routers.Get(8), Ipv6RawSocketFactory::GetTypeId());
    relay4Sock->Bind(relay4Addr);
    relay4Sock->SetAttribute("Protocol", UintegerValue(153));
    relay4Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> relay5Sock;
    Address relay5Addr(Inet6SocketAddress(interfaces[29].GetAddress(0, 1)));
    relay5Sock = Socket::CreateSocket(routers.Get(10), Ipv6RawSocketFactory::GetTypeId());
    relay5Sock->Bind(relay5Addr);
    relay5Sock->SetAttribute("Protocol", UintegerValue(153));
    relay5Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> relay6Sock;
    Address relay6Addr(Inet6SocketAddress(interfaces[31].GetAddress(0, 1)));
    relay6Sock = Socket::CreateSocket(routers.Get(11), Ipv6RawSocketFactory::GetTypeId());
    relay6Sock->Bind(relay6Addr);
    relay6Sock->SetAttribute("Protocol", UintegerValue(153));
    relay6Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> relay7Sock;
    Address relay7Addr(Inet6SocketAddress(interfaces[34].GetAddress(0, 1)));
    relay7Sock = Socket::CreateSocket(routers.Get(12), Ipv6RawSocketFactory::GetTypeId());
    relay7Sock->Bind(relay7Addr);
    relay7Sock->SetAttribute("Protocol", UintegerValue(153));
    relay7Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> relay8Sock;
    Address relay8Addr(Inet6SocketAddress(interfaces[36].GetAddress(0, 1)));
    relay8Sock = Socket::CreateSocket(routers.Get(13), Ipv6RawSocketFactory::GetTypeId());
    relay8Sock->Bind(relay8Addr);
    relay8Sock->SetAttribute("Protocol", UintegerValue(153));
    relay8Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> relay9Sock;
    Address relay9Addr(Inet6SocketAddress(interfaces[40].GetAddress(0, 1)));
    relay9Sock = Socket::CreateSocket(routers.Get(16), Ipv6RawSocketFactory::GetTypeId());
    relay9Sock->Bind(relay9Addr);
    relay9Sock->SetAttribute("Protocol", UintegerValue(153));
    relay9Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> relay10Sock;
    Address relay10Addr(Inet6SocketAddress(interfaces[46].GetAddress(0, 1)));
    relay10Sock = Socket::CreateSocket(routers.Get(22), Ipv6RawSocketFactory::GetTypeId());
    relay10Sock->Bind(relay10Addr);
    relay10Sock->SetAttribute("Protocol", UintegerValue(153));
    relay10Sock->SetRecvCallback(MakeCallback(&transit));

    Ptr<Socket> reorder1Sock;
    Address reorder1Addr(Inet6SocketAddress(interfaces[2].GetAddress(0, 1)));
    reorder1Sock = Socket::CreateSocket(edges.Get(2), Ipv6RawSocketFactory::GetTypeId());
    reorder1Sock->Bind(reorder1Addr);
    reorder1Sock->SetAttribute("Protocol", UintegerValue(153));
    if (mode == "MYREORDER")
    {
        reorder1Sock->SetRecvCallback(MakeCallback(&recvwithMYRORDER));
    }
    if (mode == "NOREORDER")
    {
        reorder1Sock->SetRecvCallback(MakeCallback(&NOREORDER));
    }
    if (mode == "CS")
    {
        reorder1Sock->SetRecvCallback(MakeCallback(&CS));
    }
    if (mode == "RCBRB")
    {
        reorder1Sock->SetRecvCallback(MakeCallback(&RCBRB));
    }

    Ptr<Socket> reorder2Sock;
    Address reorder2Addr(Inet6SocketAddress(interfaces[3].GetAddress(0, 1)));
    reorder2Sock = Socket::CreateSocket(edges.Get(3), Ipv6RawSocketFactory::GetTypeId());
    reorder2Sock->Bind(reorder2Addr);
    reorder2Sock->SetAttribute("Protocol", UintegerValue(153));
    if (mode == "MYREORDER")
    {
        reorder2Sock->SetRecvCallback(MakeCallback(&recvwithMYRORDER));
    }
    if (mode == "NOREORDER")
    {
        reorder2Sock->SetRecvCallback(MakeCallback(&NOREORDER));
    }
    if (mode == "CS")
    {
        reorder2Sock->SetRecvCallback(MakeCallback(&CS));
    }
    if (mode == "RCBRB")
    {
        reorder2Sock->SetRecvCallback(MakeCallback(&RCBRB));
    }

    Ptr<Socket> recv1Sock;
    Address recv1Addr(Inet6SocketAddress(interfaces[2].GetAddress(1, 1)));
    recv1Sock = Socket::CreateSocket(hosts.Get(2), Ipv6RawSocketFactory::GetTypeId());
    recv1Sock->Bind(recv1Addr);
    recv1Sock->SetAttribute("Protocol", UintegerValue(153));
    recv1Sock->SetRecvCallback(MakeCallback(&recv2));

    Ptr<Socket> recv2Sock;
    Address recv2Addr(Inet6SocketAddress(interfaces[3].GetAddress(1, 1)));
    recv2Sock = Socket::CreateSocket(hosts.Get(3), Ipv6RawSocketFactory::GetTypeId());
    recv2Sock->Bind(recv2Addr);
    recv2Sock->SetAttribute("Protocol", UintegerValue(153));
    recv2Sock->SetRecvCallback(MakeCallback(&recv2));

    Ptr<Socket> resend1Sock;
    Address resend1Addr(Inet6SocketAddress(interfaces[0].GetAddress(0, 1)));
    resend1Sock = Socket::CreateSocket(hosts.Get(0), Ipv6RawSocketFactory::GetTypeId());
    resend1Sock->Bind(resend1Addr);
    resend1Sock->SetAttribute("Protocol", UintegerValue(153));
    resend1Sock->SetRecvCallback(MakeCallback(&resend));

    Ptr<Socket> resend2Sock;
    Address resend2Addr(Inet6SocketAddress(interfaces[1].GetAddress(0, 1)));
    resend2Sock = Socket::CreateSocket(hosts.Get(1), Ipv6RawSocketFactory::GetTypeId());
    resend2Sock->Bind(resend2Addr);
    resend2Sock->SetAttribute("Protocol", UintegerValue(153));
    resend2Sock->SetRecvCallback(MakeCallback(&resend));

    std::string chunkEid = "20000000000000000000";
    Ptr<UniformRandomVariable> pathSelector = CreateObject<UniformRandomVariable>();
    std::vector<Address> allRelayAddresses = {relay1Addr, relay2Addr, relay3Addr, relay4Addr, relay5Addr, relay6Addr, relay7Addr, relay8Addr, relay9Addr, relay10Addr};

    double alpha_start = 0.706;
    double alpha_mid = 0.604;
    double alpha_end = 0.706;
    double beta = 92.0;
    int num_cycle = 30;
    int half_cycle = num_cycle / 2;
    double minFlowSize = 64 * 128;
    double maxFlowSize = 1024 * 1024 * 2;
    int minPktSize = 64;
    int maxPktSize = 1024;

    double baseTime = 5.0;
    double intervalBetweenPackets = 0.1;

    for (int i = 0; i < 50; i++)
    {
        int h = 0;
        if (i % 2 == 0)
        {
            h = 0;
        }
        else
        {
            h = 1;
        }

        for (int j = 0; j < 10; j++)
        {

            double alpha;
            if (j < half_cycle)
            {
                alpha = alpha_start + (alpha_mid - alpha_start) * (j / (double)half_cycle);
            }
            else
            {
                alpha = alpha_mid + (alpha_end - alpha_mid) * ((j - half_cycle) / (double)(num_cycle - half_cycle));
            }

            chunkEid = stringAddition(chunkEid, "1");
            double randValue = pathSelector->GetValue(0, 1);
            int numPaths;
            if (randValue <= 0)
            {
                numPaths = 1;
            }
            else if (randValue <= 0.1)
            {
                numPaths = 2;
            }
            else if (randValue <= 0.2)
            {
                numPaths = 3;
            }
            else if (randValue <= 0.3)
            {
                numPaths = 4;
            }
            else if (randValue <= 0.5)
            {
                numPaths = 5;
            }
            else if (randValue <= 0.6)
            {
                numPaths = 6;
            }
            else if (randValue <= 0.7)
            {
                numPaths = 7;
            }
            else if (randValue <= 0.8)
            {
                numPaths = 8;
            }
            else if (randValue <= 0.9)
            {
                numPaths = 9;
            }
            else
            {
                numPaths = 10;
            }
            std::random_shuffle(allRelayAddresses.begin(), allRelayAddresses.end());
            double flowSize = generateParetoFlowSize(alpha, minFlowSize, maxFlowSize);
            int pktSize = calculatePktSize(flowSize, minPktSize, maxPktSize);
            int pktCount = static_cast<int>(flowSize / pktSize);
            if (pktCount == 0)
                pktCount = 1;

            cout << pktCount << "," << pktSize << "\n";

            for (int k = 0; k < numPaths; k++)
            {

                Ptr<Socket> sendsock = Socket::CreateSocket(hosts.Get(h), Ipv6RawSocketFactory::GetTypeId());
                sendsock->SetAttribute("Protocol", UintegerValue(153));
                Ptr<MySend> send = CreateObject<MySend>();
                send->Setup(sendsock, allRelayAddresses[k], pktSize, pktCount / numPaths, DataRate("200Mbps"), chunkEid, recveids[h], k * (pktCount / numPaths) + 1);
                hosts.Get(0)->AddApplication(send);
                double startTime = baseTime + i * intervalBetweenPackets;
                send->SetStartTime(Seconds(startTime));
                send->SetStopTime(Seconds(startTime + 1.0));
            }
        }
    }

/*
  std::vector<FlowInfo> flowData = ReadFlowDataFromFile("flow_data.txt");

    if (flowData.empty()) {
        std::cerr << "No data read from file." << std::endl;
        return 1;
    }

    NodeContainer hosts;
    hosts.Create(2);
    Ipv6InterfaceContainer allRelayAddresses;

    double baseTime = 0.0;
    double intervalBetweenPackets = 0.1;

    for (size_t i = 0; i < flowData.size(); ++i) {
        const FlowInfo &flow = flowData[i];

        std::random_shuffle(allRelayAddresses.Begin(), allRelayAddresses.End());
        for (int k = 0; k < flow.numPaths; ++k) {
            Ptr<Socket> sendsock = Socket::CreateSocket(hosts.Get(0), Ipv6RawSocketFactory::GetTypeId());
            sendsock->SetAttribute("Protocol", UintegerValue(153));
            Ptr<MySend> send = CreateObject<MySend>();
            send->Setup(sendsock, allRelayAddresses[k], flow.pktSize, flow.pktCount / flow.numPaths, 
                        DataRate("200Mbps"), flow.eid, "some_recveid", k * (flow.pktCount / flow.numPaths) + 1);
            hosts.Get(0)->AddApplication(send);

            double startTime = baseTime + i * intervalBetweenPackets;
            send->SetStartTime(Seconds(startTime));
            send->SetStopTime(Seconds(startTime + 1.0));
        }
    }

    Simulator::Run();
    Simulator::Destroy();

*/

    if (pcap)
    {
        // p2p[11].EnablePcap("test5", devices[11].Get(0), true);
        //    p2p[0].EnablePcap("test1", devices[0].Get(0), true);
        //      p2p[0].EnablePcap("test1", devices[0].Get(1), true);
        //    p2p[9].Enaboutfile2lePcap("test9", devices[9].Get(0), true);
        //     p2p[11].EnablePcap("test2", devices[3].Get(0), true);
    }
    NS_LOG_INFO("Start Simulation.");
    Simulator::Stop(Seconds(120));
    Simulator::Run();
    NS_LOG_INFO("Finish Simulation.");
    outfile3.close();
    outfile4.close();
}
