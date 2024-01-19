#include <fstream>
#include <iostream>
#include <limits>
#include <cassert>

#include <QApplication>
#include <QCheckBox>
#include <QColor>
#include <QDir>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QPainter>
#include <QPushButton>
#include <QSpinBox>
#include <QTextStream>

#include "ui.hpp"


using namespace CDT;

static QString fileDir = QString("/Users/chenxinyu/Desktop/mesh generation/code/cdt2d/test/inputs");
static int MIN_VTX_IDX = 2;
static int MIN_EDG_IDX = -1;

CDTWidget::CDTWidget(QWidget* parent)
     :QWidget(parent)
     ,m_isDisplayIndices(true)
     ,m_delaunayRefine(true)
     ,m_translation(0., 0.)
     ,m_scale(1.0)
     ,m_sceneBox(Box2D(0,0,0,0))
     ,curInsertVtxIdx(MIN_VTX_IDX)
     ,curInsertEdgeIdx(MIN_EDG_IDX)
{
    setAutoFillBackground(false);
}

QSize CDTWidget::sizeHint() const
{
    return QSize(9001, 9001); // over 9000!
}

void CDTWidget::buildCDT(QListWidgetItem* item)
{
    QDir dir = QDir(fileDir, tr("*.txt"));
    const QString fileName = dir.filePath(item->text());
    readData(fileName);
    initTransform();
    updateCDT();
}

void CDTWidget::displayIndices(int isDisplayIndices)
{
    m_isDisplayIndices = (isDisplayIndices != 0);
    update();
}

void CDTWidget::delaunayRefine(int isDelaunayRefine)
{
   m_delaunayRefine = (isDelaunayRefine != 0);
   updateCDT();
}

void CDTWidget::prtScn()
{
    QFile file("cdt_screenshot.png");
    file.open(QIODevice::WriteOnly);
    QPixmap pixmap(rect().size());
    pixmap.fill(Qt::transparent);
    paint_(&pixmap);
    pixmap.save(&file, "PNG");
}

void CDTWidget::saveToOff()
{
    std::ofstream fout("out.off");
    fout.precision(std::numeric_limits<PrecisionType>::digits10 + 1);
    if(!fout.is_open())
        throw std::runtime_error("Save can't open file for writing OFF");
    fout << "OFF\n";
    const std::vector<Vertex>& vertices = m_cdt.GetVertices();
    const std::vector<Triangle>& triangles = m_cdt.GetTriangles();
    // Write number of vertices and faces
    fout << vertices.size() << ' ' << triangles.size() << " 0\n";
    // Write vertices
    const Box2D box = calBox2D(m_points);
    for (const Vertex& vtx : vertices)
    {
        fout << vtx.X() << ' ' << vtx.Y() << "\n";
    }
    // Write faces
    for (const Triangle& tri : triangles)
    {
       fout << "3 " << tri.V1() << ' ' << tri.V2() << ' ' << tri.V3() << "\n";
    }
    fout.close();
}

void CDTWidget::insertAllVerticesEdges()
{
   updateCDT();
   curInsertVtxIdx = MIN_VTX_IDX;
   curInsertEdgeIdx = MIN_EDG_IDX;
}

void CDTWidget::insertVertex()
{
   curInsertVtxIdx++;
   updateInsertVertex();
}

void CDTWidget::undoInsertVertex()
{
   curInsertVtxIdx--;
   updateInsertVertex();
}

void CDTWidget::updateInsertVertex()
{
   m_cdt = Triangulation();
   curInsertVtxIdx = std::min(curInsertVtxIdx, int(m_points.size()) - 1);
   curInsertVtxIdx = std::max(curInsertVtxIdx, MIN_VTX_IDX);
   std::vector<Vertex> insertVertices(m_points.begin(), m_points.begin() + curInsertVtxIdx + 1);
   m_cdt.InsertVertices(insertVertices);
   m_cdt.Perform(false);
   update();
}

void CDTWidget::insertEdge()
{
   if (m_edges.empty())
      return;
   curInsertEdgeIdx++;
   updateInsertEdge();
}

void CDTWidget::undoInsertEdge()
{
   if (m_edges.empty())
      return;
   curInsertEdgeIdx--;
   updateInsertEdge();
}

void CDTWidget::updateInsertEdge()
{
   m_cdt = Triangulation();
   m_cdt.InsertVertices(m_points);
   curInsertEdgeIdx = std::min(curInsertEdgeIdx, int(m_edges.size()) - 1);
   curInsertEdgeIdx = std::max(curInsertEdgeIdx, MIN_EDG_IDX);
   if (curInsertEdgeIdx >= 0)
   {
      std::vector<Edge> insertEdges(m_edges.begin(), m_edges.begin() + curInsertEdgeIdx + 1);
      m_cdt.InsertEdges(insertEdges);
   }
   m_cdt.Perform(false);
   update();
}

void CDTWidget::readData(const QString& fileName)
{
    QFile file(fileName);
    if(!file.open(QFile::ReadOnly))
    {
        QMessageBox::warning(
            this, tr("CDT"), tr("Could not open file ") + fileName);
    }
    QTextStream inStream(&file);
    std::size_t nPts, nEdges;
    inStream >> nPts >> nEdges;
    m_points.clear();
    for(std::size_t i = 0; i < nPts; ++i)
    {
        PrecisionType x1, y1;
        inStream >> x1 >> y1;
        m_points.emplace_back(x1, y1);
    }
    assert(m_points.size() > MIN_VTX_IDX);
    m_edges.clear();
    for(std::size_t i = 0; i < nEdges; ++i)
    {
        IdxType v1, v2;
        inStream >> v1 >> v2;
        m_edges.push_back(Edge(v1, v2));
    }
    inStream.skipWhiteSpace();
}

void CDTWidget::updateCDT()
{
    m_cdt = Triangulation();
    if(!m_points.empty())
    {
        //const CDT::DuplicatesInfo duplInfo = CDT::RemoveDuplicates(m_points);
        //if(!duplInfo.duplicates.empty())
        //{
        //    QMessageBox errBox;
        //    errBox.setText(QStringLiteral("Triangulation has duplicates"));
        //    errBox.exec();
        //}
        m_cdt.InsertVertices(m_points);
        if(!m_edges.empty())
        {
            //CDT::RemapEdges(m_edges, duplInfo.mapping);
            m_cdt.InsertEdges(m_edges);
        }
        m_cdt.Perform(m_delaunayRefine);
    }
    update();
}

void CDTWidget::paintEvent(QPaintEvent*)
{
    paint_(this);
}

QPointF CDTWidget::sceneToScreen(const Vertex& xy) const
{
    const QPointF screenCenter(width() / 2.0, height() / 2.0);
    return QPointF(m_scale * xy.X(), -m_scale * xy.Y()) + screenCenter +
           m_translation;
}

QPointF CDTWidget::screenToScene(const QPointF& xy) const
{
    const QPointF screenCenter(width() / 2.0, height() / 2.0);
    QPointF out = (xy - m_translation - screenCenter) / m_scale;
    out.setY(-out.y());
    return out;
}

double CDTWidget::calculateScale(const int w, const int h) const
{
    const Vertex sceneSize(
        m_sceneBox.XMax() - m_sceneBox.XMin(),
        m_sceneBox.YMax() - m_sceneBox.YMin()
    );
    const double sceneRatio = sceneSize.X() / sceneSize.Y();
    const double screenRatio = static_cast<double>(w) / h;
    double scale =
        (sceneRatio > screenRatio) ? w / sceneSize.X() : h / sceneSize.Y();
    return scale * 0.95;
}

void CDTWidget::initTransform()
{
    m_sceneBox = calBox2D(m_points);
    QPointF sceneCenter(
        (m_sceneBox.XMin() + m_sceneBox.XMax()) / PrecisionType(2),
        (m_sceneBox.YMin() + m_sceneBox.YMax()) / PrecisionType(2)
    );
    m_scale = calculateScale(width(), height());
    m_translation =
        QPointF(-m_scale * sceneCenter.x(), m_scale * sceneCenter.y());
}

void CDTWidget::paint_(QPaintDevice* pd)
{
    QPainter p(pd);
    p.setRenderHints(QPainter::Antialiasing);
    p.setBrush(QBrush(Qt::white));
    const std::vector<Vertex>& vertices = m_cdt.GetVertices();
    const std::vector<Triangle>& triangles = m_cdt.GetTriangles();
    const EdgeUSet& constrainedEdges = m_cdt.GetConstrainedEdges();
    if(vertices.empty())
        return;
    QPen pen;
    pen.setCapStyle(Qt::RoundCap);
    // Draw triangles
    pen.setWidthF(2.0);
    pen.setColor(QColor(150, 150, 150));
    p.setPen(pen);
    for (const Triangle& tri : triangles)
    {
        const Vertex& v1 = vertices[tri.V1()];
        const Vertex& v2 = vertices[tri.V2()];
        const Vertex& v3 = vertices[tri.V3()];
        std::vector<QPointF> pts(3);
        pts[0] = sceneToScreen(v1);
        pts[1] = sceneToScreen(v2);
        pts[2] = sceneToScreen(v3);
        const QPointF c(
           (pts[0].x() + pts[1].x() + pts[2].x()) / 3.f,
           (pts[0].y() + pts[1].y() + pts[2].y()) / 3.f);
        p.drawPolygon(pts.data(), static_cast<int>(pts.size()));
    }
    // display triangles indices
    if(m_isDisplayIndices)
    {
        int i = 0;
        for (const Triangle& tri : triangles)
        {
            const Vertex& v1 = vertices[tri.V1()];
            const Vertex& v2 = vertices[tri.V2()];
            const Vertex& v3 = vertices[tri.V3()];
            std::vector<QPointF> pts(3);
            pts[0] = sceneToScreen(v1);
            pts[1] = sceneToScreen(v2);
            pts[2] = sceneToScreen(v3);
            const QPointF c(
               (pts[0].x() + pts[1].x() + pts[2].x()) / 3.f,
               (pts[0].y() + pts[1].y() + pts[2].y()) / 3.f);
            p.drawText(c, QString::number(i));
            i++;
        }
    }
    // constraint edges
    pen.setColor(QColor(50, 50, 50));
    p.setPen(pen);
    for (const Edge& edge : constrainedEdges)
    {
        const Vertex& v1 = vertices[edge.V1()];
        const Vertex& v2 = vertices[edge.V2()];
        p.drawLine(sceneToScreen(v1), sceneToScreen(v2));
    }
    // draw points
    pen.setColor(QColor(3, 102, 214));
    pen.setWidthF(7.0);
    p.setPen(pen);
    for (const Vertex& vtx : vertices)
    {
        const QPointF pos = sceneToScreen(vtx);
        p.drawPoint(pos);
    }
    // draw points indices
    if(m_isDisplayIndices)
    {
        pen.setColor(QColor(100, 100, 255));
        p.setPen(pen);
        int i = 0;
        for (const Vertex & vtx : vertices)
        {
            const QPointF pos = sceneToScreen(vtx);
            p.drawText(pos, QString::number(i));
            i++;
        }
    }
}

void CDTWidget::mousePressEvent(QMouseEvent* event)
{
    m_prevMousePos = event->pos();
    qApp->setOverrideCursor(Qt::ClosedHandCursor);
    setMouseTracking(true);
}

void CDTWidget::mouseMoveEvent(QMouseEvent* event)
{
    m_translation += (event->pos() - m_prevMousePos);
    m_prevMousePos = event->pos();
    update();
}

void CDTWidget::mouseReleaseEvent(QMouseEvent*)
{
    qApp->restoreOverrideCursor();
    setMouseTracking(false);
}

void CDTWidget::wheelEvent(QWheelEvent* event)
{
    const double newScale =
        m_scale * std::max(0.3, (1. + event->delta() * 8e-4));
    if(m_scale == newScale)
    {
        return;
    }
    const QPointF cursor = event->posF();
    const QPointF scenePt = screenToScene(cursor);
    const QPointF screenCenter = QPointF(width(), height()) / 2.0;
    m_translation = cursor - newScale * QPointF(scenePt.x(), -scenePt.y()) -
                    screenCenter;
    m_scale = newScale;
    update();
}

void CDTWidget::resizeEvent(QResizeEvent* e)
{
    const double scaleRatio =
        calculateScale(width(), height()) /
        calculateScale(e->oldSize().width(), e->oldSize().height());
    m_scale *= scaleRatio;
    m_translation *= scaleRatio;
    update();
}

MainWindow::MainWindow(QWidget* parent)
    : QWidget(parent)
{
    m_cdtWidget = new CDTWidget();
    m_cdtWidget->setMinimumSize(QSize(1, 1));
    // Right pane
    QListWidget* filesList = new QListWidget();
    connect(
        filesList,
        SIGNAL(itemDoubleClicked(QListWidgetItem*)),
        m_cdtWidget,
        SLOT(buildCDT(QListWidgetItem*)));
    QCheckBox* displayIndices =
        new QCheckBox(QStringLiteral("Display point/triangle indices"));
    connect(
        displayIndices,
        SIGNAL(stateChanged(int)),
        m_cdtWidget,
        SLOT(displayIndices(int)));
    m_cdtWidget->displayIndices(0);
    displayIndices->setChecked(true);
    QCheckBox* delaunayRefine =
       new QCheckBox(QStringLiteral("Delaunay refinement"));
    connect(
       delaunayRefine,
       SIGNAL(stateChanged(int)),
       m_cdtWidget,
       SLOT(delaunayRefine(int)));
    m_cdtWidget->delaunayRefine(0);
    delaunayRefine->setChecked(true);
    QPushButton* insertVerticesEdgesBtn = new QPushButton(tr("Insert all vertices and edges"));
    connect(insertVerticesEdgesBtn, SIGNAL(clicked()), m_cdtWidget, SLOT(insertAllVerticesEdges()));
    insertVerticesEdgesBtn->setMinimumHeight(50);
    QPushButton* insertVtxBtn = new QPushButton(tr("Insert a vertex"));
    connect(insertVtxBtn, SIGNAL(clicked()), m_cdtWidget, SLOT(insertVertex()));
    insertVtxBtn->setMinimumHeight(50);
    QPushButton* undoInsertVtxBtn = new QPushButton(tr("Undo insert a vertex"));
    connect(undoInsertVtxBtn, SIGNAL(clicked()), m_cdtWidget, SLOT(undoInsertVertex()));
    undoInsertVtxBtn->setMinimumHeight(50);
    QPushButton* insertEdgeBtn = new QPushButton(tr("Insert a edge"));
    connect(insertEdgeBtn, SIGNAL(clicked()), m_cdtWidget, SLOT(insertEdge()));
    insertEdgeBtn->setMinimumHeight(50);
    QPushButton* undoInsertEdgeBtn = new QPushButton(tr("Undo insert a edge"));
    connect(undoInsertEdgeBtn, SIGNAL(clicked()), m_cdtWidget, SLOT(undoInsertEdge()));
    undoInsertEdgeBtn->setMinimumHeight(50);
    QGridLayout* rightLayout = new QGridLayout;
    int cntr = 0;
    rightLayout->addWidget(filesList, cntr++, 0);
    rightLayout->addWidget(displayIndices, cntr++, 0);
    rightLayout->addWidget(delaunayRefine, cntr++, 0);
    rightLayout->addWidget(insertVerticesEdgesBtn, cntr++, 0);
    rightLayout->addWidget(insertVtxBtn, cntr++, 0);
    rightLayout->addWidget(undoInsertVtxBtn, cntr++, 0);
    rightLayout->addWidget(insertEdgeBtn, cntr++, 0);
    rightLayout->addWidget(undoInsertEdgeBtn, cntr++, 0);
    // Center
    QHBoxLayout* centralLayout = new QHBoxLayout;
    centralLayout->addWidget(m_cdtWidget);
    centralLayout->addLayout(rightLayout);
    setLayout(centralLayout);
    setWindowTitle(tr("CDT Visualizer"));
    // Read files list
    QDir dir = QDir(fileDir, tr("*.txt"));
    QFileInfoList list = dir.entryInfoList();
    filesList->clear();
    QFileInfoList::iterator it;
    for(it = list.begin(); it != list.end(); ++it)
        filesList->addItem(it->fileName());
    filesList->setCurrentRow(0);
}
