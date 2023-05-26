#pragma once

#include <fstream>
#include <iostream>
#include <limits>
#include <cassert>

#include <QListWidget>
#include <QPaintEvent>

#include "cdt2d.hpp"


using namespace CDT;


struct Box2D
{
   explicit Box2D(PrecisionType xMax, PrecisionType yMax, PrecisionType xMin, PrecisionType yMin)
   {
      max.first = xMax;
      max.second = yMax;
      min.first = xMin;
      min.second = yMin;
   }
   PrecisionType XMax() const { return max.first; };
   PrecisionType YMax() const { return max.second; };
   PrecisionType XMin() const { return min.first; };
   PrecisionType YMin() const { return min.second; };
private:
   std::pair<PrecisionType, PrecisionType> max;
   std::pair<PrecisionType, PrecisionType> min;
};


static Box2D calBox2D(const std::vector<Vertex>& pts)
{
   PrecisionType xMax = pts[0].X();
   PrecisionType xMin = pts[0].X();
   PrecisionType yMin = pts[0].Y();
   PrecisionType yMax = pts[0].Y();
   for (const Vertex& pt: pts)
   {
      xMin = std::min(pt.X(), xMin);
      xMax = std::max(pt.X(), xMax);
      yMin = std::min(pt.Y(), yMin);
      yMax = std::max(pt.Y(), yMax);
   }
   return Box2D(xMax, yMax, xMin, yMin);
}

class CDTWidget : public QWidget
{
    Q_OBJECT

public:
    explicit CDTWidget(QWidget* parent = Q_NULLPTR);
    QSize sizeHint() const;

public slots:
    void buildCDT(QListWidgetItem* item);
    void displayIndices(int isDisplayIndices);
    void delaunayRefine(int isDelaunayRefine);
    void prtScn();
    void saveToOff();
    void insertAllVerticesEdges();
    void insertVertex();
    void undoInsertVertex();
    void updateInsertVertex();
    void insertEdge();
    void undoInsertEdge();
    void updateInsertEdge();

private:
    void readData(const QString& fileName);
    void updateCDT();

protected:
    void paintEvent(QPaintEvent*);

private:
    QPointF sceneToScreen(const Vertex& xy) const;
    QPointF screenToScene(const QPointF& xy) const;
    double calculateScale(const int w, const int h) const;
    void initTransform();
    void paint_(QPaintDevice* pd);
    void mousePressEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent*);
    void wheelEvent(QWheelEvent* event);
    void resizeEvent(QResizeEvent* e);

private:
    Triangulation m_cdt;
    std::vector<Vertex> m_points;
    std::vector<Edge> m_edges;
    bool m_isRemoveOuter;
    bool m_isDisplayIndices;
    bool m_delaunayRefine;
    int curInsertVtxIdx;
    int curInsertEdgeIdx;

    QPointF m_prevMousePos;
    QPointF m_translation;
    double m_scale;
    Box2D m_sceneBox;
};

class MainWindow : public QWidget
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = Q_NULLPTR);

private:
    CDTWidget* m_cdtWidget;
};

