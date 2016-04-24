#include "storage.h"

Storage::Storage(const QString& name)
    :m_name(name)
{

}

QString Storage::name() const
{
    return m_name;
}
