#ifndef LBS_GROUP
#define LBS_GROUP

//################################################################### Class def
namespace lbs
{
  /**Object holding a grouping.*/
  class LBSGroup
  {
  public:
    int id;

  public:
    LBSGroup() : id(-1) {}
    explicit LBSGroup(int in_id) : id(in_id) {}
  };
}

#endif