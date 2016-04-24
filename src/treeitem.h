#ifndef TREEITEM_H
#define TREEITEM_H

#include <QObject>
#include <QVector>
#include <QVariant>

class TreeItem
{
public:
    explicit TreeItem(const QVector<QVariant> &data, TreeItem *parentItem = 0);
    ~TreeItem();

    void appendChild(TreeItem *child);

    TreeItem *child(int row);
    int childCount() const;
    int columnCount() const;
    int childNumber() const;//function is used to determine the index of the child in its parent's list of children
    QVariant data(int column) const;
    int row() const;
    TreeItem *parentItem();

    bool insertChildren(int position, int count, int columns);
    bool insertColumns(int position, int columns);
    bool removeChildren(int position, int count);
    bool removeColumns(int position, int columns);
    bool setData(int column, const QVariant &value);

private:
    QVector<TreeItem*> m_childItems;
    QVector<QVariant> m_itemData;
    TreeItem *m_parentItem;
};
#endif // TREEITEM_H
