package muset.util;

import java.util.HashMap;
import java.util.List;
import java.util.Map;



public class SupportFunctions
{

  public static <T> Map<T,Integer> invert(List<T> list)
  {
    Map<T,Integer> result = new HashMap<T,Integer>();
    for (int i =0; i < list.size(); i++)
      result.put(list.get(i),i);
    return result;
  }

}
