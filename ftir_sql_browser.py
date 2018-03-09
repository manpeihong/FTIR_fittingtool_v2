
try:
	import MySQLdb
except:
	print( "Need to install mysql plugin, run: pip install mysqlclient")
	exit()
import sqlite3
import configparser
import os
import sys
import datetime
from PyQt5 import QtCore, QtGui, uic, QtWidgets, QtSql
from PyQt5.QtWidgets import QAbstractScrollArea
from PyQt5.QtCore import QSettings, QThread
import numpy as np

os.chdir(os.path.dirname(os.path.realpath(sys.argv[0])))  # Change the working directory to current directory.
qtDesignerFile = "ftir_sql_browser.ui" # GUI layout file.
Ui_MainWindow, QtBaseClass = uic.loadUiType( qtDesignerFile )


class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):
	# Create signals here
	signal_to_emit = QtCore.pyqtSignal(str)
	def __init__(self):
		QtWidgets.QMainWindow.__init__(self)
		Ui_MainWindow.__init__(self)
		self.setupUi(self)

		config = configparser.ConfigParser()
		config.read('configuration.ini')

		self.Connect_To_SQL( config )
		self.sqlUser_lineEdit.setText( config["SQL_Server"]["default_user"] );

		self.refresh_commandLinkButton.clicked.connect( self.Reinitialize_Tree_Table )
		self.sqlUser_lineEdit.returnPressed.connect( self.Reinitialize_Tree_Table )
		self.filter_lineEdit.textEdited.connect( self.Reinitialize_Tree_Table )
		self.finished_pushButton.clicked.connect( self.close )

		self.meta_data = {}
		self.background_x_data, self.background_y_data = [], []
		self.x_data, self.y_data = [], []

		#settings = QSettings( "configuration.ini", QSettings.IniFormat, self )
		#while( not self.Initialize_DB_Connection( settings.value( "SQL_Server/database_type" ),
		#							  settings.value( "SQL_Server/host_location" ),
		#							  settings.value( "SQL_Server/database_name" ),
		#							  settings.value( "SQL_Server/username" ),
		#							  settings.value( "SQL_Server/password" ) ) ):
		#	print( "Error with SQL Connection" )
		#	delay_seconds = 10;
		#	print( "Trying again in {} seconds".format( delay_seconds ) )
		#	QThread.sleep( delay_seconds )
		#self.sqlUser_lineEdit.setText( settings.value( "SQL_Server/default_user" ) );

		self.Initialize_Tree_Table()

		## Can only connect functions after previous things are initialized
		#self.my_pushButton.pressed.connect( lambda : self.signal_to_emit.emit( "You pushed the button" ) )
		#self.signal_to_emit.connect( self.Pop_Up_Message )

	def Connect_To_SQL( self, configuration_file ):
		if configuration_file['SQL_Server']['database_type'] == "QSQLITE":
			try:
				self.sql_db = sqlite3.connect( configuration_file['SQL_Server']['database_name'] )
			except sqlite3.Error as e:
				error = QtWidgets.QMessageBox()
				error.setIcon( QtWidgets.QMessageBox.Critical )
				error.setText( e )
				error.setWindowTitle( "Unable to connect to SQL Database" )
				error.exec_()
				print( e )
				self.quit()
				return
		elif configuration_file['SQL_Server']['database_type'] == "QMYSQL":
			try:
				self.sql_db = MySQLdb.connect(host=configuration_file['SQL_Server']['host_location'],db=configuration_file['SQL_Server']['database_name'],
									user=configuration_file['SQL_Server']['username'],passwd=configuration_file['SQL_Server']['password'])
				self.sql_db.ping( True ) # Maintain connection to avoid timing out
			except:
				error = QtWidgets.QMessageBox()
				error.setIcon( QtWidgets.QMessageBox.Critical )
				#error.setText( e )
				error.setWindowTitle( "Unable to connect to SQL Database" )
				error.exec_()
				#print( e )
				self.quit()
				return


	#def Initialize_DB_Connection( self, database_type, host_location, database_name, username, password ):
	#	self.sql_db = QtSql.QSqlDatabase.addDatabase( database_type );
	#	if database_type == "QMYSQL":
	#		self.sql_db.setConnectOptions( "MYSQL_OPT_RECONNECT=1" );

	#	if not host_location:
	#		self.sql_db.setHostName( host_location )
	#	self.sql_db.setDatabaseName( database_name )
	#	if not username:
	#		self.sql_db.setUserName( username )
	#	if not password:
	#		self.sql_db.setPassword( password )
	#	sql_worked = self.sql_db.open()
	#	if sql_worked:
	#		print( "Connected to {} Database {}".format( host_location, database_name ) )
	#	else:
	#		problem = self.sql_db.lastError()
	#		print( "Error with SQL Connection: ", problem )

	#	return sql_worked


	def Initialize_Tree_Table( self ):
		self.Reinitialize_Tree_Table()
		#		connect( ui.treeWidget, &QTreeWidget::itemDoubleClicked, [this]( QTreeWidgetItem* tree_item, int column )
		#{
		#	vector<const QTreeWidgetItem*> things_to_graph = this->Get_Bottom_Children_Elements_Under( tree_item );

		#	for( const QTreeWidgetItem* x : things_to_graph )
		#	{
		#		Graph_Tree_Node( x );
		#	}
		#} );

		# setup policy and connect slot for context menu popup:
		self.treeWidget.setContextMenuPolicy( QtCore.Qt.CustomContextMenu );
		self.treeWidget.customContextMenuRequested.connect( self.treeContextMenuRequest );

	def Reinitialize_Tree_Table( self ):
		self.treeWidget.clear()
		self.treeWidget.setSizeAdjustPolicy( QAbstractScrollArea.AdjustToContents )

		user = self.sqlUser_lineEdit.text()
		header_titles = [ "Sample Name", "Date", "Temperature (K)", "Bias (V)", "Time of Day", "measurement_id" ]
		self.treeWidget.setHeaderLabels( header_titles )
		self.treeWidget.hideColumn( len( header_titles ) - 1 )
		what_to_collect = [ "sample_name", "date(time)", "temperature_in_k", "bias_in_v", "time(time)", "measurement_id" ]

		first_filter = "sample_name LIKE \"%" + self.filter_lineEdit.text() + "%\""
		user_filter = "user = \"" + user + "\""

		self.Recursive_Tree_Table_Build( what_to_collect, self.treeWidget.invisibleRootItem(), 0, [ first_filter, user_filter ] )
		for i in range( len(header_titles) ):
			self.treeWidget.resizeColumnToContents( i )

	def Recursive_Tree_Table_Build( self, what_to_collect, parent_tree, current_collectable_i, filters ):
		if( current_collectable_i == len( what_to_collect ) ):
			return

		query = self.sql_db.cursor()
		querey_string = "SELECT DISTINCT {} FROM ftir_measurements".format( what_to_collect[ current_collectable_i ] )

		if filters:
			requirements = " AND ".join( filters );
			querey_string += " WHERE {}".format( requirements );

		try:
			test = query.execute( querey_string )
			selected_rows = query.fetchall()
		except MySQLdb.Error as e:
			print( "Error pulling data from ftir_measurments:%d:%s" % (e.args[0], e.args[1]) )
			return

		numberOfRows = len( selected_rows )

		for row in selected_rows:
			current_value = row[ 0 ] # Only grabbing from one column
			if isinstance( current_value, datetime.timedelta ):
				hours, remainder = divmod(current_value.total_seconds(), 3600)
				minutes, seconds = divmod(remainder, 60) 
				current_value = '%02d:%02d:%02d' % (hours, minutes, seconds)

			new_filters = filters.copy()
			if current_value is None:
				new_filters.append( "{} IS NULL".format( what_to_collect[ current_collectable_i ] ) )
			else:
				new_filters.append( "{} = \"{}\"".format( what_to_collect[ current_collectable_i ], str( current_value ) ) )
			new_tree_branch = parent_tree; # Only add a new breakout for the first one, and if more than 1 child
			if numberOfRows > 1 or current_collectable_i == 0:
				new_tree_branch = QtWidgets.QTreeWidgetItem( parent_tree )

			new_tree_branch.setText( current_collectable_i, str( current_value ) );
			self.Recursive_Tree_Table_Build( what_to_collect, new_tree_branch, current_collectable_i + 1, new_filters )

	def Get_Bottom_Children_Elements_Under( self, tree_item ):
		number_of_children = tree_item.childCount()
		if number_of_children == 0:
			return [tree_item]

		lowest_level_children = []
		for i in range( number_of_children ):
			i_lowest_level_children = self.Get_Bottom_Children_Elements_Under( tree_item.child( i ) )
			lowest_level_children += i_lowest_level_children

		return lowest_level_children

	def treeContextMenuRequest( self, pos ):
		menu = QtWidgets.QMenu( self )
		menu.setAttribute( QtCore.Qt.WA_DeleteOnClose );
		selected = self.treeWidget.selectedItems()
		actually_clicked = self.treeWidget.itemAt( pos )
		if len( selected ) == 1 and len( self.Get_Bottom_Children_Elements_Under( selected[0] ) ) == 1: # context menu on legend requested
			menu.addAction( "Set As Data", lambda : self.Set_Sample_Data( selected[ 0 ].text( selected[ 0 ].columnCount() - 1 ) ) )
			menu.addAction( "Set As Background", lambda : self.Set_Background( selected[ 0 ].text( selected[ 0 ].columnCount() - 1 ) ) )

		menu.popup( self.treeWidget.mapToGlobal( pos ) )

	def Grab_SQL_Data_From_Measurement_ID( self, measurement_id ):
		query = self.sql_db.cursor()
		if not query.execute( "SELECT wavenumber,intensity FROM raw_ftir_data WHERE measurement_id = \"{}\"".format(measurement_id) ):
			#qDebug() << "Error pulling data from ftir_measurments: "
			#	<< query.lastError();
			return [[],[]]
		selected_rows = query.fetchall()
		return zip(*selected_rows) # Unfolds into x and y components

	def Grab_SQL_Metadata_From_Measurement_ID( self, measurement_id ):
		query = self.sql_db.cursor()
		what_to_collect = [ "sample_name", "date(time)", "temperature_in_k", "bias_in_v", "time(time)" ]
		if not query.execute( "SELECT DISTINCT {} FROM ftir_measurements WHERE measurement_id = \"{}\"".format( ','.join( what_to_collect ), measurement_id ) ):
			#qDebug() << "Error pulling data from ftir_measurments: "
			#	<< query.lastError();
			return {}
		selected_rows = query.fetchall()
		if selected_rows:
			return { k:v for k,v in zip(what_to_collect, selected_rows[0]) } # Builds dictionary with results with what_to_collect strings as the keys
		else:
		   return {}

	def Set_Background( self, measurement_id ):
		background_meta_data = self.Grab_SQL_Metadata_From_Measurement_ID( measurement_id )
		self.selectedBackground_label.setText( ' '.join( [str(v) for v in background_meta_data.values()] ) )
		self.background_x_data, self.background_y_data = self.Grab_SQL_Data_From_Measurement_ID( measurement_id )
	def Set_Sample_Data( self, measurement_id ):
		self.meta_data = self.Grab_SQL_Metadata_From_Measurement_ID( measurement_id )
		self.selectedSample_label.setText( ' '.join( [str(v) for v in self.meta_data.values()] ) )
		self.x_data, self.y_data = self.Grab_SQL_Data_From_Measurement_ID( measurement_id )

	def Get_Meta_Data( self ):
		return self.meta_data
	def Get_Background_Data( self ):
		return self.background_x_data, self.background_y_data
	def Get_Sample_Data( self ):
		return self.x_data, self.y_data

def Find_Nearest_Value_In( look_in_x_y, x_value ):
	idx = ( np.abs( np.array(look_in_x_y[0]) - x_value ) ).argmin()
	return look_in_x_y[1][idx]

def Remove_Background( sample, background ):
	final_data = [sample[0],[]]
	for x, y in zip( *sample ):
		final_data[1].append( y / Find_Nearest_Value_In( background, x ) )

	return final_data

def Get_Data():
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	#sys.exit(app.exec_())
	app.exec_()

	meta_data = window.Get_Meta_Data()
	background = window.Get_Background_Data()
	sample = window.Get_Sample_Data()

	return meta_data, Remove_Background( sample, background )
