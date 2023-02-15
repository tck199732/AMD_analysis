import pathlib
import pickle
import sqlite3
import pandas as pd

class DatabaseSQL:
    def __init__(self, fname, columns, df_name='dataframe', primary_key='id', overwrite=False):
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
        overwrite : bool
            If True, the table will be created even if it already exists.
        """

        self.fname = fname
        self.df_name = df_name
        self.primary_key = primary_key

        if not df_name in columns.keys():
            columns[df_name] = 'BLOB'

        if primary_key in columns.keys():
            ptype = columns[primary_key]
        else:
            ptype = 'INTEGER'
        ptype = f'{ptype} PRIMARY KEY'

        self.columns = columns
        self.connection = sqlite3.connect(self.fname)
        self.cursor = self.connection.cursor()

        # create table with columns
        execute_str = ',\n'.join([
            f'{name} {type}' for name, type in columns.items() if not name == primary_key
        ])

        execute_str = f'{primary_key} {ptype},\n' + execute_str
        print(f'CREATE TABLE IF NOT EXISTS dataframes\n({execute_str})')
        if overwrite:
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

    def _exists(self, **kwargs):
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

    def add_dataframe(self, dataframe, forced_update=True, **kwargs):
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

        if self._exists(**kwargs):
            if not forced_update:
                return
            query = f'UPDATE dataframes SET {self.df_name}=? ' + ' '.join([
                f'AND {name}=?' for name in kwargs.keys()
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

        if self._exists(**kwargs):
            query = f'SELECT {self.df_name} FROM dataframes WHERE ' + ' '.join([
                f'{name}=? AND' for name in kwargs.keys()
            ])
            values = tuple(kwargs.values())
            dataframe_bytes = self.cursor.execute(query, values).fetchone()[0]
            dataframe = pickle.loads(dataframe_bytes)
        else:
            dataframe = None

        return dataframe