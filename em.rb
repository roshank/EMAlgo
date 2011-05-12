require 'logger'

#Intialization is done by picking random numbers in K slices
#Termination done on stabalizing ALL means to within .001 OR hitting iteration maximum

class EM
  class Matrix
    def initialize(inputs, means, init_type)
      @logger = Logger.new(STDOUT)
      @logger.level = Logger::WARN
      @x = inputs
      @mu = means
      if (init_type == 1)
        distributed_initialization()
      else 
        random_initialization
      end
      @z = Array.new(@x.length) { Array.new(@mu.length, 0.0) }
    end
    
    # Disribute starting points randomly but in evenly spaced slices
    def distributed_initialization()
      slices = @x.max / @mu.length
      (0..@mu.length-1).each do |j|
        @mu[j] = slices*j + rand()*(slices)
      end
    end
    
    # distribute starting points randomly within range of values 0...max(input)
    def random_initialization()
      (0..@mu.length-1).each do |j|
        @mu[j] = rand() * @x.max
      end
    end
    
    
    # do the last printout of each Zij
    def final_print
      prob = Array.new(@x.length, 0.0)
      (0...@mu.length).each do |j|
        (0...@x.length).each do |i|
          if (i >= 25)
            break
          end
          # calc p(x1 | z11=1), p(x2 | z21=1), etc
          @z[i][j] = Math.exp( 0 - ((@x[i] - @mu[j])**2)/2)
          p "p(x#{i+1} | z#{i+1}#{j+1}) = #{@z[i][j]}"
          prob[i] += @z[i][j]
        end
      end    
    end
    
    # E step
    def expectation()
      prob = Array.new(@mu.length, 0.0)
      (0...@mu.length).each do |j|
        (0...@x.length).each do |i|
          
          # calc p(x1 | z11=1), p(x2 | z21=1), etc
          @z[i][j] = Math.exp( 0 - ((@x[i] - @mu[j])**2)/2)
          @logger.debug("z[#{i+1}][#{j+1}] = Math.exp( 0 - ((#{@x[i]} - #{@mu[j]})**2)/2) =  #{@z[i][j]}")
        end
        
        # Finished a complete row of inputs (zi1, for ex)
        # so calculate the total prob of this row
        total_prob = 0.0
        (0..@x.length-1).each do |i|
          total_prob += @z[i][j]
          prob[j] += @z[i][j]
        end
        
        # to_print =""
        # Divide the probabilty of each element with the total prob
        (0..@x.length-1).each do |i|
          @z[i][j] = @z[i][j] / total_prob
          #to_print = "#{to_print} @z[#{i+1}][#{j+1}]#{@z[i][j].to_s},"
        end
        
      end
      
      # Multiply probability of each X having occured to determine probability of this set of data from this Mu 
      total = 1.0
      (0..prob.length-1).each do |x|
        total = total * prob[x]
        total = Math.log(total.abs)
        p "Chance data sampled from Mu#{x+1} is #{total}"
      end
      
    end
    
    
    # M step
    def maximization()
      mu_diff = 0.0
      (0..@mu.length-1).each do |j|
        total_prob = 0.0
        value = 0.0
        (0..@x.length-1).each do |i|
          total_prob += @z[i][j]
          value += @z[i][j] * @x[i]
        end
        mu_diff = [mu_diff, (@mu[j] - value / total_prob).abs].max
        @mu[j] = value / total_prob
        p "mu[#{j+1}]=#{@mu[j]},"
      end
      return mu_diff
    end
  end  
  
  def initialize(inputs, num_clusters, init_type)
    @inputs = inputs
    @means = Array.new(num_clusters, 0)
    @em = Matrix.new(@inputs, @means, init_type)
  end
  
  def run(n)
    (0..n).each do |i|
      p "Iteration #{i}"
      @em.expectation()
      mu_diff = @em.maximization()
      if (mu_diff < 0.001)
        p "Terminating because MEAN values have stabalized after #{i} iterations"
        break
      end
      if (i == n)
        p "Terminating because max iterations reached"
        break
      end
    end
    @em.final_print()
    
    #BIC = 2 * ln L(x | θ-hat) - r ln n
    # L ( x | 0-hat ) - num_clusters * ln inputs.length
    # What is theta-hat??
    
    bic = 0 
    
    #bic = 2 * LL - Klusters * Log(n);
    
    p "BIC Score: #{bic}"
  end
  
end

def in_from_file(path)
  arr = []
  arr = File.open(path).read.split.map{|s| s.to_f }
  return arr
end


if __FILE__ == $0
  arr = [9, 10, 11, 20, 21, 22, 46, 49, 55, 57]
  (1..5).each do |k|
    p "Executing EM with K = #{k}"
    num_clusters = k
    x = EM.new(arr, num_clusters, 1)
    x.run(100)
    p "Completed run of EM with K = #{k}"
  end
  
  arr = in_from_file('in.txt')
  (1..5).each do |k|
    p "Executing EM with K = #{k}"
    num_clusters = k
    x = EM.new(arr, num_clusters, 1)
    x.run(100)
    p "Completed run of EM with K = #{k}"
  end  
end
