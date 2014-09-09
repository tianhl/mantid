#ifndef MULTIDATASETFIT_H_
#define MULTIDATASETFIT_H_

#include "MantidQtAPI/UserSubWindow.h"
#include "MantidAPI/AlgorithmObserver.h"
#include "MantidQtAPI/WorkspaceObserver.h"
#include "ui_MultiDatasetFit.h"
#include "ui_AddWorkspace.h"

#include <QMap>
#include <vector>

class QwtPlot;
class QTableWidget;
class QComboBox;
class QPushButton;

namespace MantidQt
{
namespace CustomInterfaces
{

class PlotController;
class DatasetPlotData;

/**
 * Class MultiDatasetFitDialog implements a dialog for setting up a multi-dataset fit
 * and displaying the results.
 */

class MultiDatasetFit: public API::UserSubWindow, public Mantid::API::AlgorithmObserver,
                          public MantidQt::API::WorkspaceObserver
{
  Q_OBJECT
public:
  /// The name of the interface as registered into the factory
  static std::string name() { return "Multi dataset fitting"; }
  // This interface's categories.
  static QString categoryInfo() { return "General"; }
  /// Constructor
  MultiDatasetFit(QWidget *parent = NULL);
  ~MultiDatasetFit();

signals:
  void dataTableUpdated();

protected:
  /// To be overridden to set the appropriate layout
  virtual void initLayout();

private slots:
  void addWorkspace();
  void workspaceSelectionChanged();
  void removeSelectedSpectra();
private:
  void addWorkspaceSpectrum(const QString &wsName, int wsIndex);
  /// The form generated by Qt Designer
  Ui::MultiDatasetFit m_uiForm;
  /// Controls the plot and plotted data.
  PlotController *m_plotController;
};

/*==========================================================================================*/
/**
  * A dialog for selecting a workspace from the ADS.
  */
class AddWorkspaceDialog: public QDialog
{
  Q_OBJECT
public:
  AddWorkspaceDialog(QWidget *parent);
  QString workspaceName() const {return m_workspaceName;} 
  std::vector<int> workspaceIndices() const {return m_wsIndices;}
private slots:
  void accept();
  void reject();
  void workspaceNameChanged(const QString&);
  void selectAllSpectra(int state);
private:
  /// Name of the selected workspace
  QString m_workspaceName;
  /// Selected workspace index
  std::vector<int> m_wsIndices;
  /// Maximum index in the selected workspace
  int m_maxIndex;
  Ui::AddWorkspace m_uiForm;
};

/*==========================================================================================*/
/**
  * A class for controlling the plot widget and the displayed data.
  */
class PlotController: public QObject
{
  Q_OBJECT
public:
  PlotController(QObject *parent, QwtPlot *plot, QTableWidget *table, QComboBox *plotSelector, QPushButton *prev, QPushButton *next);
  ~PlotController();
  void clear();
private slots:
  void tableUpdated();
  void prevPlot();
  void nextPlot();
  void plotDataSet(int);
private:
  /// The plot widget
  QwtPlot *m_plot;
  /// The workspace table
  QTableWidget *m_table;
  QComboBox *m_plotSelector;
  QPushButton *m_prevPlot;
  QPushButton *m_nextPlot;
  QMap<int,boost::shared_ptr<DatasetPlotData>> m_plotData;
  int m_currentIndex;
};

} // CustomInterfaces
} // MantidQt

#endif /*MULTIDATASETFITDIALOG_H_*/
