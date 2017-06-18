Node = Struct.new(:x, :y, :next_node)

class ListReader
  attr_accessor :first_node

  def initialize
    @first_node = nil
    load_file("files/list.txt")
  end

  def load_file(filename)
    prev_node = nil
    begin
      File.open(filename, 'r') do |f|
        while line = f.gets
          values = line.split(',')
          node = Node.new(values[0].to_i, values[1].to_f, nil)
          if prev_node
            prev_node.next_node = node
            prev_node = node
          else
            self.first_node = node
            prev_node = node
          end
        end
      end
    rescue Exception => e
      puts "There was an error opening the file, make sure is located in the same folder as the program and that is called list.txt"
    end
  end
end

class MathCalculator
  def self.sum(list)
    list.inject(0.0){|sum, e| sum + e}
  end

  def self.mean(list)
    sum(list)/list.size
  end

  def self.standard_deviation(list)
    Math.sqrt(_sum_of_difference(list)/(list.size-1))
  end

  def self._sum_of_difference(list)
    avg = mean(list)
    list.inject(0.0){|sum, e| sum += (e - avg)**2}
  end
end

class TDistribution
  E = 0.0000001

  def self.calculate_p(w, x, num_seg, dof)
    (w/3.0)*(f(0, dof) + odd_sum(num_seg, w, dof) + even_sum(num_seg, w, dof) + f(x, dof))
  end

  def self.even_sum(num_seg, w, dof)
    (2..num_seg-2).step(2).inject(0){|result, element| result + (2*f(element*w, dof))}
  end

  def self.f(x, dof)
    temp = (dof+1)/2.0
    (gamma(temp)/(((dof*Math::PI)**(0.5))*gamma(dof/2.0)))*((1+((x**2)/dof.to_f))**(-1*temp))
  end

  def self.gamma(x)
    return 1 if x == 1
    return Math.sqrt(Math::PI) if x == 0.5
    return (x - 1)*gamma(x - 1)
  end

  def self.odd_sum(num_seg, w, dof)
    (1..num_seg-1).step(2).inject(0){|result, element| result + (4*f(element*w, dof))}
  end

  def self.p_value(t, dof)
    num_seg = 20
    w = t/num_seg
    current_p = calculate_p(w, t, num_seg/2, dof)
    new_p = calculate_p(w, t, num_seg, dof)
    while((current_p - new_p).abs > E)
      num_seg *= 2
      w = t/num_seg
      current_p = new_p
      new_p = calculate_p(w, t, num_seg, dof)
    end
    new_p
  end

  def self.t_value(desired_p, dof)
    t = 1.0
    delta = 0.5
    p_value = p_value(t, dof)
    current_error = p_value - desired_p
    delta_sign = current_error > 0 ? 1 : -1
    while current_error.abs > 0.0000001
      t = delta_sign > 0 ? t - delta : t + delta
      p_value = p_value(t, dof)
      current_error = p_value - desired_p
      current_error_sign = current_error > 0 ? 1 : -1
      if current_error_sign != delta_sign
        delta = delta / 2
        delta_sign = current_error_sign
      end
    end
    return t
  end
end

class Probe
  attr_reader :list_of_x, :list_of_y
  attr_reader :sum_of_x, :sum_of_y, :sum_of_xy, :sum_of_x2, :sum_of_y2
  attr_reader :mean_of_x, :mean_of_y, :n, :estimated_proxy

  def initialize(list_of_x, list_of_y, estimated_proxy = 60.17)
    @list_of_x = list_of_x
    @list_of_y = list_of_y
    @n = list_of_x.size
    @estimated_proxy = estimated_proxy
    @sum_of_x = MathCalculator.sum(list_of_x)
    @sum_of_y = MathCalculator.sum(list_of_y)
    @sum_of_xy = MathCalculator.sum(_list_of_xy)
    @sum_of_x2 = MathCalculator.sum(_list_of_squares(list_of_x))
    @sum_of_y2 = MathCalculator.sum(_list_of_squares(list_of_y))
    @mean_of_x = MathCalculator.mean(list_of_x)
    @mean_of_y = MathCalculator.mean(list_of_y)
  end

  def _list_of_squares(list)
    list.inject([]){|sum, e| sum << e**2}
  end

  def _list_of_xy
    result = []
    list_of_x.each_with_index do |e, i|
      result << e * list_of_y[i]
    end
    result
  end

  def beta_1
    return 0 if n == 0
    (sum_of_xy - (n * mean_of_x * mean_of_y)) / (sum_of_x2 - (n * (mean_of_x**2)))
  end

  def beta_0
    mean_of_y - (beta_1 * mean_of_x)
  end

  def r
    return 0 if n == 0
    ((n * sum_of_xy) - (sum_of_x * sum_of_y)) / Math.sqrt(((n * sum_of_x2) - (sum_of_x**2)) * ((n * sum_of_y2) - (sum_of_y**2)))
  end

  def r_square
    r * r
  end

  def _standar_deviation_term
    result = []
    list_of_x.each_with_index do |e, i|
      result << (list_of_y[i] - beta_0 - (beta_1 * e))**2
    end
    Math.sqrt((1.0/(n - 2)) * MathCalculator.sum(result))
  end

  def _range_helper
    Math.sqrt(1 + (1.0/n) + ((estimated_proxy - MathCalculator.mean(list_of_x))/MathCalculator._sum_of_difference(list_of_x)))
  end

  def range
    TDistribution.t_value(0.35, n - 2) * _standar_deviation_term * _range_helper
  end

  def upi
    projected_a_and_m_size + range
  end

  def lpi
    projected_a_and_m_size - range
  end

  def significance
    t = (r.abs * Math.sqrt(n-2))/Math.sqrt(1 - r_square)
    p_value = TDistribution.p_value(t, n-2)
    return 1 - (2 * p_value)
  end

  def projected_a_and_m_size
    beta_0 + beta_1 * estimated_proxy
  end
end

class ResultPrinter
  attr_accessor :list_reader

  def initialize(list_reader)
    @list_reader = list_reader
  end

  def print_result
    puts "|#{'Parameter'.center(15)}|#{'Result'.center(13)}|"
    node = list_reader.first_node
    list_of_x = []
    list_of_y = []
    while(node)
      list_of_x << node.x
      list_of_y << node.y
      node = node.next_node
    end
    probe = Probe.new(list_of_x, list_of_y)
    puts "|#{'r:'.center(15)}| #{probe.r.round(9)} |"
    puts "|#{'r2:'.center(15)}| #{probe.r_square.round(9).to_s.center(11)} |"
    puts "|#{'Significance:'.center(15)}| #{probe.significance.round(11)} |"
    puts "|#{'B0:'.center(15)}| #{probe.beta_0.round(9)} |"
    puts "|#{'B1:'.center(15)}| #{probe.beta_1.round(9)} |"
    puts "|#{'P:'.center(15)}| #{probe.projected_a_and_m_size.round(8)} |"
    puts "|#{'Range:'.center(15)}| #{probe.range} |"
    puts "|#{'UPI:'.center(15)}| #{probe.upi.round(8)} |"
    puts "|#{'LPI:'.center(15)}| #{probe.lpi.round(7)} |"
  end
end

#MainCode
list = ListReader.new
ResultPrinter.new(list).print_result if list.first_node