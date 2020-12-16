import os
import pandas
import numpy as np
import datetime
from datetime import datetime
from os.path import join as pjoin
from data_loader import EventsDataset


class EmailDataset(EventsDataset):

    def __init__(self, split, data_dir="./Email/"):
        super(EmailDataset, self).__init__()

        filepath = data_dir+"sub-1d/"
        if split == 'train':
            time_start = 0
            # time_end = datetime(1971, 3, 31, tzinfo=self.TZ).toordinal()
            #time_end = datetime(1971, 2, 28, tzinfo=self.TZ).toordinal()
            time_end = datetime(1970, 5, 31, tzinfo=self.TZ).toordinal()
        elif split == 'test':
            # time_start = datetime(1971, 4, 1, tzinfo=self.TZ).toordinal()
            # time_end = datetime(1971, 6, 10, tzinfo=self.TZ).toordinal()
            #time_start = datetime(1971, 3, 1, tzinfo=self.TZ).toordinal()
            #time_end = datetime(1971, 6, 10, tzinfo=self.TZ).toordinal()
            time_start = datetime(1970, 6, 1, tzinfo=self.TZ).toordinal()
            time_end = datetime(1970, 8, 1, tzinfo=self.TZ).toordinal()
        else:
            raise ValueError('invalid split', split)

        self.FIRST_DATE = datetime(1970, 1, 1, tzinfo=self.TZ)

        '''self.TEST_TIMESLOTS = [datetime(1970, 3, 1, tzinfo=self.TZ),
                               datetime(1970, 3, 20, tzinfo=self.TZ),
                               datetime(1971, 4, 10, tzinfo=self.TZ),
                               datetime(1971, 5, 1, tzinfo=self.TZ),
                               datetime(1971, 5, 20, tzinfo=self.TZ),
                               datetime(1971, 6, 10, tzinfo=self.TZ)]'''
        self.TEST_TIMESLOTS = [datetime(1970, 6, 10, tzinfo=self.TZ),
                               datetime(1970, 6, 20, tzinfo=self.TZ),
                               datetime(1970, 6, 30, tzinfo=self.TZ),
                               datetime(1970, 7, 10, tzinfo=self.TZ),
                               datetime(1970, 7, 20, tzinfo=self.TZ),
                               datetime(1970, 8, 1, tzinfo=self.TZ)]

        data_file = pjoin(data_dir, 'email-Eu-core-temporal-Dept3_sub.csv')
        csv = pandas.read_csv(data_file, encoding='utf-8')
        self.data = {}
        to_date1 = lambda s: datetime.strptime(s, '%Y-%m-%d')
        to_date2 = lambda s: datetime.strptime(s, '%Y-%m-%d %H:%M:%S')
        self.user_columns = list(filter(lambda c: c.find('id') >= 0, list(csv.keys())))
        assert len(self.user_columns) == 2, (list(csv.keys()), self.user_columns)
        self.time_column = list(filter(lambda c: c.find('date') >= 0, list(csv.keys())))
        assert len(self.time_column) == 1, (list(csv.keys()), self.time_column)
        self.time_column = self.time_column[0]

        # #20 node missing, total number of nodes = 89, but have index range=0-89, so built 90*90 matrix
        # self.N_nodes = len(np.unique(csv[self.user_columns[0]].tolist() + csv[self.user_columns[1]].tolist()))
        self.N_nodes = 90

        # print(np.unique(csv[self.user_columns[0]].tolist()))
        # print(np.unique(csv[self.user_columns[1]].tolist()))
        print("Number of users:%d" % self.N_nodes)


        values = []
        # for tt in range(len(csv[self.time_column])):
        #     print(csv[self.time_column][tt])
        #     values.append(csv[self.time_column][tt])
        # print(values)
        for column in list(csv.keys()):
            values = csv[column].tolist()
            for fn in [int, to_date2]:
                try:
                    values = list(map(fn, values))
                    break
                except Exception as e:
                    continue
            self.data[column] = values
        n_rows = len(self.data[self.time_column])
        # print(type(self.data[self.time_column][0]))
        time_stamp_days = np.array([d.toordinal() for d in self.data[self.time_column]], dtype=np.int)
        conditions = [~np.isnan(self.data[self.user_columns[0]]),
                      ~np.isnan(self.data[self.user_columns[1]]),
                      np.array(self.data[self.user_columns[0]]) != np.array(self.data[self.user_columns[1]]),
                      time_stamp_days >= time_start,
                      time_stamp_days <= time_end]
        valid_ids = np.ones(n_rows, dtype=np.bool)
        # print("Before+++"+str(np.sum(valid_ids)))
        for cond in conditions:
            valid_ids = valid_ids & cond
        # print("After+++"+str(np.sum(valid_ids)))
        
        self.valid_ids = np.where(valid_ids)[0]
        time_stamps_secs = [self.data[self.time_column][i].timestamp() for i in self.valid_ids]
        self.valid_ids = self.valid_ids[np.argsort(time_stamps_secs)]
        for column in list(csv.keys()):
            values = csv[column].tolist()
            self.data[column] = [values[i] for i in self.valid_ids]
            
            # print(len(self.data[key]))
        


        self.A_initial = np.random.randint(0, 2, size=(self.N_nodes, self.N_nodes))
        self.A_last = np.random.randint(0, 2, size=(self.N_nodes, self.N_nodes))

        # print('\nA_initial', np.sum(self.A_initial))
        # print('A_last', np.sum(self.A_last), '\n')

        self.n_events = len(self.valid_ids)
        self.all_events = []
        for i in range(self.n_events):
            user_id1 = self.data[self.user_columns[0]][i]
            user_id2 = self.data[self.user_columns[1]][i]
            event_time = to_date2(self.data[self.time_column][i])
            self.all_events.append((user_id1, user_id2, 'communication event', event_time))

        self.event_types = ['communication event']

        self.all_events = sorted(self.all_events, key=lambda t: t[3].timestamp())
        print('\n%s' % split.upper())
        print('%d events between %d users loaded' % (len(self.all_events), self.N_nodes))
        # print(self.all_events[3])


        self.event_types_num = {'association event': 0}
        k = 1  # k >= 1 for communication events
        for t in self.event_types:
            self.event_types_num[t] = k
            k += 1
        self.n_events = len(self.all_events)

        initial_embeddings = np.zeros((90,36))
        for m in range(36):
            MotifMatrix = filepath+"%d.txt" % m
            # MotifMatrix = "./Email/sub-1h/%d.txt" % m
            motifMatrix = np.loadtxt(MotifMatrix, delimiter=" ")
            for i in range(90):
                initial_embeddings[i][m] = np.sum(motifMatrix[i])
        np.savetxt(filepath+"nodeMotifcount.txt", initial_embeddings, fmt='%.0f')
        # np.savetxt("./Email/sub-1h/nodeMotifcount.txt", initial_embeddings, fmt='%.0f')
        m = np.mean(initial_embeddings)
        mx = np.max(initial_embeddings)
        mn = np.min(initial_embeddings)
        for m in range(36):
            for i in range(90):
                initial_embeddings[i][m] = (float(initial_embeddings[i][m])-m) / (mx - mn)

        self.initial_embeddings = initial_embeddings
        np.savetxt(filepath+"nodeMotifcount_initialembedding.txt", initial_embeddings, fmt='%.3f')
        # np.savetxt("./Email/sub-1h/nodeMotifcount_initialembedding.txt", initial_embeddings, fmt='%.3f')

    def get_Adjacency(self, multirelations=False):
        if multirelations:
            print('warning: this dataset has only one relation type, so multirelations are ignored')
        return self.A_initial, ['association event'], self.A_last
