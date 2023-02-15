import pathlib
import pickle
import sqlite3
import pandas as pd

class DatabaseSQL:
    def __init__(self, fname, columns=None, df_name='dataframe', primary_key='id', mode='a'):
        """ Initialize a database with a table and columns.
        Parameters
        ----------
        fname : str
            The name of the database file.
        columns : dict
            A dictionary of column names and their types.
        df_name : str
            The name of the column that will store dataframes.
        primary_key : str
            The name of the primary key column.
        mode : str
            The mode to open the database in. Can be 'a' for append or 'w' for write or 'r' for read.
        """

        
        self.fname = fname
        self.mode = mode
        self.df_name = df_name
        self.primary_key = primary_key
        self.columns = columns

        self.connection = sqlite3.connect(self.fname)
        self.cursor = self.connection.cursor()

        if mode == 'w' or mode == 'a':
            if not df_name in self.columns.keys():
                self.columns[df_name] = 'BLOB'

            if primary_key in self.columns.keys():
                ptype = self.columns[primary_key]
            else:
                ptype = 'INTEGER'
            ptype = f'{ptype} PRIMARY KEY'

            # create table with columns
            execute_str = ',\n'.join([
                f'{name} {type}' for name, type in columns.items() if not name == primary_key
            ])

            execute_str = f'{primary_key} {ptype},\n' + execute_str
            if mode == 'w':
                self.cursor.execute('DROP TABLE IF EXISTS dataframes')
            
            self.cursor.execute(f'CREATE TABLE IF NOT EXISTS dataframes\n({execute_str})')

        
        
    def _valid_columns(self, **kwargs):
        """ Check if kwargs are valid columns.
        Parameters
        ----------
        kwargs : dict
            The values of the columns to check.
        """
        for key in kwargs.keys():
            if not key in self.columns.keys():
                raise ValueError(f'Invalid column name: {key}')

    def exists(self, **kwargs):
        """ Check if a row exists in the database.
        Parameters
        ----------
        kwargs : dict
            The values of the columns to check.
        """

        self._valid_columns(**kwargs)
            
        query_str = ' AND '.join([
            f'{name}=?' for name in kwargs.keys()
        ])
        query = f'SELECT {self.primary_key} FROM dataframes WHERE {query_str}'
        values =  tuple(kwargs.values())
        existing_row = self.cursor.execute(query, values).fetchone()

        return existing_row is not None

    def add_dataframe(self, dataframe, forced_update=False, **kwargs):
        """ Add a dataframe to the database.
        Parameters
        ----------
        dataframe : pandas.DataFrame
            The dataframe to add to the database.
        forced_update : bool
            If True, the dataframe will be updated if it already exists.
        kwargs : dict
            The values of the columns to add to the database.
        """
        # convert dataframe to bytes
        dataframe_bytes = pickle.dumps(dataframe)

        if self.exists(**kwargs):
            if not forced_update:
                return
            query = f'UPDATE dataframes SET {self.df_name}=? ' + ' AND '.join([
                f'{name}=?' for name in kwargs.keys()
            ])
            values = (dataframe_bytes, ) + tuple(kwargs.values())
        else:
            keys = ', '.join([
                f'{name}' for name in kwargs.keys()
            ])
            query = f'INSERT INTO dataframes ({self.df_name}, {keys}) VALUES (?, ' + ', '.join(['?'] * len(kwargs)) + ')'
            values = (dataframe_bytes, ) + tuple(kwargs.values())

        self.cursor.execute(query, values)
        self.connection.commit()

    def get_dataframe(self, **kwargs):
        """ Get a dataframe from the database.
        Parameters
        ----------
        kwargs : dict
            The values of the columns to check.
        """
        try:
            query = f'SELECT {self.df_name} FROM dataframes WHERE ' + ' AND '.join([
                f'{name}=?' for name in kwargs.keys()
            ])
            values = tuple(kwargs.values())
            dataframe_bytes = self.cursor.execute(query, values).fetchone()[0]
            dataframe = pickle.loads(dataframe_bytes)
        except:
            dataframe = None

        return dataframe