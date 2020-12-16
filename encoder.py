import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
import torch.nn.init as init
import torch.nn.functional as F


# Social
tor_M = np.zeros((36, 84, 84))
initial_embeddings = np.zeros((84,36))
filepath = "./SocialEvolution/10min/"
for m in range(36):
    MotifMatrix = "./SocialEvolution/10min/%d.txt" % m
    motifMatrix = np.loadtxt(MotifMatrix, delimiter=" ")
    for i in range(84):
        initial_embeddings[i][m] = np.sum(motifMatrix[i])
np.savetxt(filepath+"nodeMotifcount.txt", initial_embeddings, fmt='%.0f')
MotifCount = "./SocialEvolution/10min/socialevol-timewindow-10min-motif-Rcounts.txt"
MC = np.loadtxt(MotifCount, delimiter=" ")

# Email
# tor_M = np.zeros((36, 90, 90))
# for m in range(36):
#     MotifMatrix = "./Email/sub-1d/%d.txt" % m
#     motifMatrix = np.loadtxt(MotifMatrix, delimiter=" ")
# MotifCount = "./Email/sub-1d/email-Eu-core-temporal-Dept3_sub-timewindow-1d-motif-Rcounts.txt"
# MC = np.loadtxt(MotifCount, delimiter=" ")

# Github
# tor_M = np.zeros((36, 284, 284))
# for m in range(36):
#     MotifMatrix = "./Github/1d/%d.txt" % m
#     motifMatrix = np.loadtxt(MotifMatrix, delimiter=" ")
# MotifCount = "./Github/1d/Github-timewindow-1d-motif-Rcounts.txt"
# MC = np.loadtxt(MotifCount, delimiter=" ")

# College
# tor_M = np.zeros((36, 1899, 1899))
# for m in range(36):
#     MotifMatrix = "./CollegeMSG/1min/%d.txt" % m
#     motifMatrix = np.loadtxt(MotifMatrix, delimiter=" ")
# MotifCount = "./CollegeMSG/1min/CollegeMsg-timewindow-1min-motif-Rcounts.txt"
# MC = np.loadtxt(MotifCount, delimiter=" ")

def normalization(Mat):
    # L=D^-0.5*(A+I)*D^-0.5
    I = np.eye(Mat.shape[0])
    Mat = Mat + I
    degree = np.array(Mat.sum(1))
    d_hat = np.diag(np.power(degree,-0.5).flatten())
    return d_hat.dot(Mat).dot(d_hat)


class LinearEncoder(nn.Module):
    def __init__(self, n_in, n_out, sym=True):
        super(LinearEncoder, self).__init__()

        self.fc = nn.Linear(n_in * 2, n_out)
        self.sym = sym
        self.init_weights()

    def init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear) or isinstance(m, nn.Bilinear):
                nn.init.xavier_normal_(m.weight.data)
                m.bias.data.fill_(0.1)

    def node2edge(self, x):
        N = x.shape[0]
        mask = x.new(1, N).fill_(1)
        node_i = torch.nonzero(mask).repeat(1, N).view(-1, 1)
        node_j = torch.nonzero(mask).repeat(N, 1).view(-1, 1)
        if self.sym:
            triu = (node_i < node_j).squeeze()  # skip loops and symmetric connections
        else:
            triu = (node_i != node_j).squeeze()  # skip loops and symmetric connections
        idx = (node_i * N + node_j)[triu].squeeze()  # linear index
        edges = torch.cat((x[node_i[triu]],
                           x[node_j[triu]]), 1).view(int(torch.sum(triu)), -1)

        return edges, idx

    def edges2matrix(self, x, idx, N):
        edges = x.new(N * N, x.shape[1]).fill_(0)
        edges[idx] = x
        edges = edges.view(N, N, -1)
        return edges

    def forward(self, inputs):
        x = inputs  # N,n_hid
        N = x.shape[0]
        x, idx = self.node2edge(x)  # Eq. 6
        x = self.fc(x)  # Eq. 7: get edge embeddings (N,N,n_hid)
        x = self.edges2matrix(x, idx, N)  # N,N,n_hid
        if self.sym:
            x = x + x.permute(1, 0, 2)
        return x



class MLP(nn.Module):
    """Two-layer fully-connected ELU net with batch norm."""

    def __init__(self, n_in, n_hid, n_out, do_prob=0., bilinear=False, bnorm=True):
        super(MLP, self).__init__()
        self.bilinear = bilinear
        self.bnorm = bnorm
        if bilinear:
            self.fc1 = nn.Bilinear(n_in, n_in, n_hid)
        else:
            self.fc1 = nn.Linear(n_in, n_hid)
        self.fc2 = nn.Linear(n_hid, n_out)
        if bnorm:
            self.bn = nn.BatchNorm1d(n_out)
        self.dropout_prob = do_prob

        self.init_weights()

    def init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear) or isinstance(m, nn.Bilinear):
                nn.init.xavier_normal(m.weight.data)
                m.bias.data.fill_(0.1)
            elif isinstance(m, nn.BatchNorm1d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()

    def batch_norm(self, inputs):
        x = self.bn(inputs)
        return x

    def forward(self, inputs):
        # Input shape: [num_sims, num_things, num_features]
        if self.bilinear:
            x = F.elu(self.fc1(inputs[0], inputs[1]))
            x = x.view(x.size(0), -1)
        else:
            x = F.elu(self.fc1(inputs))
        x = F.dropout(x, self.dropout_prob, training=self.training)
        x = F.elu(self.fc2(x))
        if self.bnorm:
            return self.batch_norm(x)
        else:
            return x


class MLPEncoder(nn.Module):
    def __init__(self, n_in, n_hid, n_out, do_prob=0., factor=True, bilinear=False, n_stages=2, bnorm=True, sym=True):
        super(MLPEncoder, self).__init__()

        self.factor = factor
        self.bilinear = bilinear
        self.n_stages = n_stages
        self.sym = sym
        if self.sym:
            raise NotImplementedError('')

        # self.mlp1 = MLP(n_in, n_hid, n_hid, do_prob, bilinear=bilinear, bnorm=bnorm)

        
        self.mlp1 = MLP(n_in, n_hid, n_hid, do_prob, bnorm=bnorm)
        self.mlp2 = MLP(n_hid * (1 if bilinear else 2), n_hid, n_hid, do_prob, bilinear=bilinear, bnorm=bnorm)
        if n_stages == 2:
            self.mlp3 = MLP(n_hid, n_hid, n_hid, do_prob)
            if self.factor:
                self.mlp4 = MLP(n_hid * (2 if bilinear else 3), n_hid, n_hid, do_prob, bilinear=bilinear, bnorm=False)
                print("Using factor graph MLP encoder.")
            else:
                self.mlp4 = MLP(n_hid * (1 if bilinear else 2), n_hid, n_hid, do_prob, bilinear=bilinear, bnorm=False)
                print("Using MLP encoder.")
        self.fc_out = nn.Linear(n_hid, n_out)


        self.TMAGCN = GcnNet(input_dim=n_hid, hid_dim=n_hid, output_dim=n_hid)
        self.TMA_num = 4
        # motif 2.0
        # self.TMAcat = nn.Linear(n_hid * (self.TMA_num+1), n_hid)
        # motif 3.0
        self.alpha = 1
        self.TMAcat = nn.Linear(n_hid * (self.TMA_num), n_hid)
        self.bili = nn.Bilinear(n_in, n_in, n_hid)
        self.init_weights()

    def init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear) or isinstance(m, nn.Bilinear):
                nn.init.xavier_normal_(m.weight.data)
                m.bias.data.fill_(0.1)

    def edge2node(self, x):
        # NOTE: Assumes that we have the same graph across all samples.
        N = x.shape[0]
        device = x.get_device() if x.is_cuda else 'cpu'
        rel_rec = (1 - torch.eye(N, device=device)).unsqueeze(2)  # to sum up all neighbors except for self loops
        incoming = (x * rel_rec).sum(dim=1)  # N,n_hid
        return incoming / N

    def node2edge(self, x):
        # NOTE: Assumes that we have the same graph across all samples.
        N = x.shape[0]
        node_i = torch.arange(N).view(N, 1).repeat(1, N).view(-1, 1)
        # [0,0,0,1,1,1,2,2,2]
        node_j = torch.arange(N).view(N, 1).repeat(N, 1).view(-1, 1)
        # [0,1,2,0,1,2,0,1,2]
        if self.sym:
            triu = (node_i < node_j).squeeze()  # skip loops and symmetric connections
            # [0,1,1,0,0,1,0,0,0]
        else:
            triu = (node_i != node_j).squeeze()  # skip loops
            # [0,1,1,1,0,1,1,1,0]
        idx = (node_i * N + node_j)[triu].squeeze()  # linear index =(1,N*N), int(torch.sum(triu) 1 and N*N-int(torch.sum(triu) 0
        if self.bilinear:
            edges = (x[node_i[triu]], x[node_j[triu]])
        else:
            edges = torch.cat((x[node_i[triu]],
                               x[node_j[triu]]), 1).view(int(torch.sum(triu)), -1)

        return edges, idx

    def edges2matrix(self, x, idx, N):
        edges = x.new(N * N, x.shape[1]).fill_(0)
        edges[idx] = x
        edges = edges.view(N, N, -1)
        return edges

    def forward(self, inputs, edges=None):
        x = inputs  # N,n_hid  prev_embeding
        N = x.shape[0]
        x = self.mlp1(x)  # f_v^1: 2-layer ELU net per node

        x, idx = self.node2edge(x)  # Eq. 6  x:(torch.sum(triu),n_hid)

        x = self.mlp2(x)  # f_e^1: get edge embeddings (N,N,n_hid)

        if self.n_stages == 2:
            x_skip = x  # edge embeddings: N*(N-1)/2, n_hid
            x = self.edges2matrix(x, idx, N)  # N,N,n_hid

            if edges is not None:
                x_skip = self.edges2matrix(x_skip, idx, N)  # N,N,n_hid

                u, v = edges[0, 0].item(), edges[0, 1].item()

                x_skip = torch.cat((x_skip[u, v].view(1, -1), x_skip[v, u].view(1, -1)), dim=0)  # 2,n_hid

                if self.sym:
                    raise NotImplementedError('')

            if self.factor:
                x = self.edge2node(x)  # N,n_hid
                x = x[[u, v], :]  # 2,n_hid
                N = 2
                x = self.mlp3(x)  # f_v^2: 2,n_hid
                x, idx = self.node2edge(x)  # N*(N-1)/2, n_hid
                if self.bilinear:
                    x = (torch.cat((x[0].view(x[0].size(0), -1), x_skip), dim=1),
                         torch.cat((x[1].view(x[1].size(0), -1), x_skip), dim=1))  # Skip connection
                else:
                    x = torch.cat((x, x_skip), dim=1)  # Skip connection
                x = self.mlp4(x)  # N*(N-1)/2, n_hid

                # print(x.size())
                x = self.edges2matrix(x, idx, N)  # N,N,n_hid
                # print(x.size())
            else:
                x = self.mlp3(x)
                x = torch.cat((x, x_skip), dim=1)  # Skip connection
                x = self.mlp4(x)
        else:
            x = self.edges2matrix(x, idx, N)  # N,N,n_hid

        # print(x.size()) #torch.Size([2, 2, 32])
        x = self.fc_out(x)  # N,N,n_out
        if self.sym:
            x = x + x.permute(1, 0, 2)

        return x, idx



class TMA_GraphConvolution(nn.Module):
    def __init__(self, input_dim, output_dim, use_bias=True):
        super(TMA_GraphConvolution, self).__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.use_bias = use_bias
        self.weight = nn.Parameter(torch.Tensor(input_dim, output_dim))
        if self.use_bias:
            self.bias = nn.Parameter(torch.Tensor(output_dim))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        init.kaiming_uniform_(self.weight)
        if self.use_bias:
            init.zeros_(self.bias)

    def forward(self, adjacency, input_feature):
        
        support = torch.mm(input_feature, self.weight)
        N = input_feature.shape[0]
        # GCN
        # output = torch.mm(adjacency, support)

        # TMA
        # use adjacency matrix
        # adjacency = adjacency.cpu().detach().numpy()
        # no adjacency matrix, only motif matrix
        adjacency = np.zeros((N, N))

        moti_adjacency = np.zeros((N, N))
        MC_array = MC.reshape(36,1)
        MC_array = np.squeeze(MC_array)
        MC_array_index = np.argsort(MC_array)
        for i in range(36):
            moti_adjacency += tor_M[MC_array_index[-i-1]]
        adjacency = adjacency + moti_adjacency

        adjacency = normalization(adjacency)
        adjacency = torch.from_numpy(adjacency).cuda().to(torch.float32)
        output = torch.mm(adjacency, support)
        

        if self.use_bias:
            output += self.bias
        # return output
        return torch.sigmoid(output)


class TMA_GcnNet(nn.Module):
    """
    定义一个包含两层GraphConvolution的模型
    """
    def __init__(self, input_dim, hid_dim, output_dim):
        super(TMA_GcnNet, self).__init__()
        self.gcn1 = TMA_GraphConvolution(input_dim, hid_dim)
        self.gcn2 = TMA_GraphConvolution(hid_dim, output_dim)
    
    def forward(self, adjacency, feature):
        h = F.relu(self.gcn1(adjacency, feature))
        logits = self.gcn2(adjacency, h)
        return logits


class GraphConvolution(nn.Module):
    def __init__(self, input_dim, output_dim, use_bias=True):
        super(GraphConvolution, self).__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.use_bias = use_bias
        self.weight = nn.Parameter(torch.Tensor(input_dim, output_dim))
        if self.use_bias:
            self.bias = nn.Parameter(torch.Tensor(output_dim))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        init.kaiming_uniform_(self.weight)
        if self.use_bias:
            init.zeros_(self.bias)

    def forward(self, adjacency, input_feature):
        
        support = torch.mm(input_feature, self.weight)
        output = torch.mm(adjacency, support)

        if self.use_bias:
            output += self.bias
        # return output
        return torch.sigmoid(output)


class GcnNet(nn.Module):
    """
    定义一个包含两层GraphConvolution的模型
    """
    def __init__(self, input_dim, hid_dim, output_dim):
        super(GcnNet, self).__init__()
        self.gcn1 = GraphConvolution(input_dim, hid_dim)
        self.gcn2 = GraphConvolution(hid_dim, output_dim)
    
    def forward(self, adjacency, feature):
        h = F.relu(self.gcn1(adjacency, feature))
        logits = self.gcn2(adjacency, h)
        return logits

class GraphAttentionLayer(nn.Module):
    """
    Simple GAT layer, similar to https://arxiv.org/abs/1710.10903
    """
    def __init__(self, in_features, out_features, dropout, alpha, concat=True):
        super(GraphAttentionLayer, self).__init__()
        self.dropout = dropout
        self.in_features = in_features
        self.out_features = out_features
        self.alpha = alpha
        self.concat = concat

        self.W = nn.Parameter(torch.empty(size=(in_features, out_features)))
        nn.init.xavier_uniform_(self.W.data, gain=1.414)
        self.a = nn.Parameter(torch.empty(size=(2*out_features, 1)))
        nn.init.xavier_uniform_(self.a.data, gain=1.414)

        self.leakyrelu = nn.LeakyReLU(self.alpha)

    def forward(self, h, adj):        
        Wh = torch.mm(h, self.W) # h.shape: (N, in_features), Wh.shape: (N, out_features)
        a_input = self._prepare_attentional_mechanism_input(Wh)
        e = self.leakyrelu(torch.matmul(a_input, self.a).squeeze(2))

        # adj for motif
        MC_array = MC.reshape(36,1)
        MC_array = np.squeeze(MC_array)
        MC_array_index = np.argsort(MC_array)
        adjacency = adj.cpu().numpy()
        for i in range(4):
            adjacency += tor_M[MC_array_index[-i-1]]
        adj = normalization(adjacency)
        adj = torch.from_numpy(adj).cuda()
        

        zero_vec = -9e15*torch.ones_like(e)
        attention = torch.where(adj > 0, e, zero_vec)
        attention = F.softmax(attention, dim=1)
        attention = F.dropout(attention, self.dropout, training=self.training)
        h_prime = torch.matmul(attention, Wh)

        if self.concat:
            return F.elu(h_prime)
        else:
            return h_prime

    def _prepare_attentional_mechanism_input(self, Wh):
        N = Wh.size()[0] # number of nodes

        
        Wh_repeated_in_chunks = Wh.repeat_interleave(N, dim=0)
        Wh_repeated_alternating = Wh.repeat(N, 1)
        

        all_combinations_matrix = torch.cat([Wh_repeated_in_chunks, Wh_repeated_alternating], dim=1)
        # all_combinations_matrix.shape == (N * N, 2 * out_features)

        return all_combinations_matrix.view(N, N, 2 * self.out_features)

    def __repr__(self):
        return self.__class__.__name__ + ' (' + str(self.in_features) + ' -> ' + str(self.out_features) + ')'


class GAT(nn.Module):
    def __init__(self, nfeat, nhid, dropout, alpha, nheads):
        super(GAT, self).__init__()
        self.dropout = dropout

        self.attentions = [GraphAttentionLayer(nfeat, nhid, dropout=dropout, alpha=alpha, concat=True) for _ in range(nheads)]
        for i, attention in enumerate(self.attentions):
            self.add_module('attention_{}'.format(i), attention)

        self.out_att = GraphAttentionLayer(nhid * nheads, nhid, dropout=dropout, alpha=alpha, concat=False)

    def forward(self, x, adj):
        x = F.dropout(x, self.dropout, training=self.training)
        x = torch.cat([att(x, adj) for att in self.attentions], dim=1)
        x = F.dropout(x, self.dropout, training=self.training)
        x = F.elu(self.out_att(x, adj))
        return F.log_softmax(x, dim=1)