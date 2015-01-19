#ifndef MANTIDQTCUSTOMINTERFACES_DATACOMPARISON_H_
#define MANTIDQTCUSTOMINTERFACES_DATACOMPARISON_H_

//----------------------
// Includes
//----------------------
#include "ui_DataComparison.h"
#include "MantidQtAPI/UserSubWindow.h"
#include "MantidAPI/MatrixWorkspace.h"

#include <QPointer>

#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_magnifier.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_zoomer.h>


namespace MantidQt
{
namespace CustomInterfaces
{
  class DataComparison : public MantidQt::API::UserSubWindow
  {
    Q_OBJECT

  public:
    /// The name of the interface as registered into the factory
    static std::string name() { return "Data Comparison"; }
    // This interface's categories.
    static QString categoryInfo() { return "General"; }

  public:
    /// Default Constructor
    DataComparison(QWidget *parent = 0);

  private slots:
    /// Add selected data to plot
    void addData();
    /// Remove selected data from plot
    void removeSelectedData();
    /// Remove all data from plot
    void removeAllData();
    /// Create a diff of the two selected workspaces
    void diffSelected();
    /// Remove the diff from the plot
    void clearDiff();
    /// Handles replotting workspace spectra
    void plotWorkspaces();
    /// Handles updating the plot
    void updatePlot();
    /// Handles updating th eplot after a spectum index is changed
    void spectrumIndexChanged();
    /// Handles creating and plotting a diff worksapce
    void plotDiffWorkspace();
    /// Toggle the pan plot tool
    void togglePan(bool enabled);
    /// Toggle the zoom plot tool
    void toggleZoom(bool enabled);
    /// Resets the zoom level to show all curves
    void resetView();

  private:
    /// Enumeration for column index
    enum Column
    {
      COLOUR,
      WORKSPACE_NAME,
      SPEC_OFFSET,
      CURRENT_SPEC
    };

    /// Initialize the layout
    virtual void initLayout();
    /// Normalises spectra offsets in table
    void normaliseSpectraOffsets();
    /// Gets an initial curve colour for a new workspace
    int getInitialColourIndex();

  private:
    // The form generated by Qt Designer
    Ui::DataComparison m_uiForm;

    // The plot object
    QwtPlot *m_plot;
    // Curves shown on plot, indexed by workspace name
    QMap<QString, boost::shared_ptr<QwtPlotCurve>> m_curves;

    // Plot zoom tool
    QwtPlotZoomer *m_zoomTool;
    // Plot pan tool
    QwtPlotPanner *m_panTool;
    // Plot magnify tool
    QwtPlotMagnifier *m_magnifyTool;

    boost::shared_ptr<QwtPlotCurve> m_diffCurve;
    // The two workspaces that are currently being diffed
    QPair<QString, QString> m_diffWorkspaceNames;

  };

}
}

#endif //MANTIDQTCUSTOMINTERFACES_DATACOMPARISON_H_
