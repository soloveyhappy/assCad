#include "mainwindow.h"
#include <QStatusBar>
#include <QLayout>
#include <QMouseEvent>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QApplication>
#include <QDebug>
#include <QVector>
#include <QMenu>
#include <QMenuBar>

#include "treemodel.h"
#include "treeitem.h"

#define EMPTYLIST QVector<QVariant>()

MainWindow::MainWindow()
{

    m_mdiArea = new QMdiArea();
    setCentralWidget(m_mdiArea);
    m_dataDoc = new QDockWidget( tr("Data"));
    m_treeView = new QTreeView;
    m_dataDoc->setWidget(m_treeView);
    addDockWidget(Qt::DockWidgetArea::LeftDockWidgetArea, m_dataDoc);
    resize(1024, 860);
    m_dataDoc->raise();

    root = new TreeItem(QVector<QVariant>()<<"Name"<<"Info");
    group1 = new TreeItem(QVector<QVariant>()<<"group1"<<"g1");
    group2 = new TreeItem(QVector<QVariant>()<<"group2"<<"g2");
    TreeItem* subgroup1 = new TreeItem(QVector<QVariant>()<<"subgroup1"<<"s1");
    TreeItem* subgroup2 = new TreeItem(QVector<QVariant>()<<"subgroup2"<<"s1");

    group2->appendChild(subgroup1);
    group2->appendChild(subgroup2);
    root->appendChild(group1);
    root->appendChild(group2);
    treeModel = new TreeModel(root);

    //m_treeView->setSelectionBehavior(QAbstractItemView::SelectionBehavior::SelectItems);
    m_treeView->setSelectionMode(QAbstractItemView::ExtendedSelection);
    m_treeView->setDragEnabled(true);
    m_treeView->setAcceptDrops(true);
    m_treeView->setDropIndicatorShown(true);
    m_treeView->setModel(treeModel);
    m_treeView->resizeColumnToContents(0);
    m_treeView->resizeColumnToContents(1);

    QAction* action = new QAction("add", menuBar());
    connect(action, SIGNAL(triggered(bool)), this, SLOT(onAction(bool)));
    menuBar()->addActions(QList<QAction*>() << action);
}

MainWindow::~MainWindow()
{

}

void MainWindow::clicked(bool)
{
    QApplication::sendEvent(this, new QMouseEvent(QEvent::MouseButtonPress,
                                                  QPointF(10,10),
                                                  Qt::MouseButton::LeftButton,
                                                  Qt::MouseButton::LeftButton,
                                                  Qt::KeyboardModifier::NoModifier ));


}
static int iter = 0;
void MainWindow::onAction(bool)
{
    TreeItem* newItem = new TreeItem(EMPTYLIST << "NewItem" <<  QString::number(++iter));
    //QModelIndex rootIndex = treeModel->index(0, 0);
    QModelIndex group1Index = treeModel->index(0, 0);
    treeModel->insertItem(newItem, group1Index);

    //group1->appendChild(newItem);
    //m_treeView->reset();
}

void MainWindow::mousePressEvent(QMouseEvent* mouseEvent)
{

}

